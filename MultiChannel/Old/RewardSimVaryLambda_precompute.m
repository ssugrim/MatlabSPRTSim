%Computes the number of samples required to find 1 good channel in a set of
%varying size with varying lambda

%for use during cmd line run in linux
%addpath '../lib'

%Multi Threading
matlabpool open;

%% Test Tuning Params
%global Parameters


%number of channels - Note if this is too small, we might not have enough
%good channels to find. 
Channels = 2:30;

%The bin we deem as the desriable one
good_bin = 1;

%penalitly coefficent;
lambda_range = (0);

%The hard limt All tests must stop at (Lenght of our vectors, should be
%longer than eta, for some allowance).
Test_Limit = 1000;

Histogram_max = 100;

%number of Trials
Trial_limit = 10;

%% SPRT params
delta = 0.125;
split = 0.5;
p_prime = 0.2;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%% left side hypothesis threhold line computation
%tuneing parameters - center and width
p_prime_ls = p_prime;

%boundaries of accept and reject region
p_0_ls = p_prime_ls - (split * delta);
p_1_ls = p_prime_ls + ((1 - split) * delta);

%generate test lines
[L_1_ls,L_0_ls] = GenTestLines(alpha,beta,p_1_ls,p_0_ls,Test_Limit);

%% right side hypothesis threhold line computation
%tuneing parameters - center and width
p_prime_rs = 1 - p_prime_ls;

%boundaries of accept and reject region
p_0_rs = p_prime_rs - (split * delta);
p_1_rs = p_prime_rs + ((1 - split) * delta); 
%generate test lines
[L_1_rs,L_0_rs] = GenTestLines(alpha,beta,p_1_rs,p_0_rs,Test_Limit);

%% Convience names
Outer_lower_boundary = L_0_ls;
Outer_upper_boundary = L_1_rs;
Inner_lower_boundary = L_1_ls;
Inner_upper_boundary = L_0_rs;

%% Storage Vars for plot generation

Mean_local_samples_to_find_1 = zeros(length(Channels),length(lambda_range));
Mean_global_samples_to_find_1 = zeros(length(Channels),length(lambda_range));
Mean_error = zeros(length(Channels),length(lambda_range));
Mean_samples_per_channel = zeros(length(Channels),length(lambda_range));
histogram = zeros(length(Channels),length(lambda_range),Histogram_max);

%% Begin Simulation

for channel_index = 1:1:length(Channels)
    %% Channel Loop
    sprintf('On Channel set %d',Channels(channel_index))
    for lambda_index = 1:1:length(lambda_range)
        %% Lambda Loop
        sprintf('On Lambda set %f',lambda_range(lambda_index))
        
        %% Trial Storage Vars
        
        %Value of the local clock when the decsion was made
        Trail_local_samples_to_find_1 = zeros(1,Trial_limit);
        
        %Value of the global clock when the decsion was made
        Trail_global_samples_to_find_1 = zeros(1,Trial_limit);
        
        %Was there an error when the decsion was made
        Trial_error = zeros(1,Trial_limit);
        
        %Mean number of sammples perchannel during this trial
        Trial_samples_per_channel = zeros(1,Trial_limit);
        
        %The histogram for this trial
        Trial_histogram = zeros(Trial_limit,Histogram_max);
        %% computed Trail Values
        %maximal cost is the distance between the outer threshold lines
        %to maximum measurement penality.
        max_cost = pdist([Test_Limit,Outer_lower_boundary(Test_Limit);Test_Limit,Outer_upper_boundary(Test_Limit)]) + (lambda_range(lambda_index) * Test_Limit);
        
        parfor trial_index = 1:Trial_limit
        %for trial_index = 1:1:Trial_limit
            %% Trial Loop
            sprintf('On Trial %d',trial_index)
            
            %% Genrerate a Random vector of appropriate length with 1 good channel in
            %pick Channels(channel_index) - 1 channels in bin2
%             p_act_part = random('unif',p_1_ls,p_0_rs,[1,Channels(channel_index)-1]);
%             
%             %choose a random isertion point
%             insert = floor(random('unif',1,Channels(channel_index)-1));
%             
%             %Insert the good_val into the insert+1 random point
%             p_act = [p_act_part(1:insert), good_val, p_act_part(insert+1:end)];
            p_act = rand(Channels(channel_index),1);
            
            %Bin each of the actual P's. 
            p_act_bin = zeros(1,Channels(channel_index));
            for k = 1:Channels(channel_index)
                if p_act(k) > p_prime_rs
                    %above the upper limit, must belong to 3
                    p_act_bin(k) = 3;
                elseif p_act(k) < p_prime_ls
                    %below the lower limit must belong to 1
                    p_act_bin(k) = 1;
                else
                    %every thing else is in 2
                    p_act_bin(k) = 2;
                end
            end
                       
            %% Trail temprorary vars we're keeping track of all the steps now
            
            %pre compute the samples list for each channel - Doingthis for
            %uniformity when we run 3 tests. 
            samples = zeros(Channels(channel_index), Test_Limit);
            
            for k = 1:Channels(channel_index)
                samples(k,:) = random('bino',1,p_act(k),[1,Test_Limit]);
            end
            
            %vector of accumulated ones - (1,1) is the starting point for
            %each channel 
            ones_sum = zeros(Channels(channel_index), Test_Limit);
            
            for k = 1:Channels(channel_index)
                ones_sum(k,:) = cumsum(samples(k,:));
            end
            
            %One clock per channel - Number of measurements each channel
            %gets (l_k) Have to start a 1 for indexing purposes. 
            local_clocks = zeros(Channels(channel_index),Test_Limit);
                    
            %vector of decision for each channel, zero if undecided, 1,2 or
            %3 for the respective bins
            decisions = zeros(1,Channels(channel_index));
            
            %Vector of costs for each channel (e_k)
            costs = zeros(Channels(channel_index),Test_Limit);
            
            %Compute the the intial cost for each channel (should all be the
            %same). we're centering on 1,1            
            for k = 1:1:Channels(channel_index)
                costs(k,1) = DistToThresh(Outer_upper_boundary,Outer_lower_boundary,ones_sum(k,1),local_clocks(k,1),Test_Limit)
            end
%TODO there may be a more efficient representation of this (perhaps just pass the coordinates)
            
            %The path of indices we sample
            path = zeros(1,Test_limit);
            
            for global_clock = 1:1:Test_Limit
                %% The actual SPRT Loop
                
                %pick a channel
                pick_index =  pick(costs(:,global_clock));
                
                %track the choices
                path(global_clock) = pick_index;
                
                %track how many measurements a channel got and when.
                local_clocks(pick_index,global_clock) = 1;
                
                %copy all the old values
                costs(:,(global_clock + 1)) = costs(:,global_clock)
                
                % Compute curret time, and add one since we're centering on (1,1)
                current_time = sum(local_clocks(pick_index,:)) + 1;
                                               
                %update the cost and %check for decisons
                if  decisions(pick_index) == 0
                    if ones_sum(pick_index,current_time) > Outer_upper_boundary()
                        %decide region 3
                        costs(pick_index,(global_clock + 1)) = max_cost;
                        decisions(k) = 3;
                    elseif ones_sum(pick_index,current_time) < Outer_lower_boundary(local_clocks(k))
                        %decide region 1
                        costs(pick_index,(global_clock + 1)) = max_cost;
                        decisions(k) = 1;
                        %we foud a good one stop the loop
                        break;
                    elseif Inner_lower_boundary(local_clocks(k)) < ones_sum(pick_index,current_time) && ones_sum(pick_index,current_time) < Inner_upper_boundary(local_clocks(k))
                        %decide region 2
                        costs(pick_index,(global_clock + 1)) = max_cost;
                        decisions(k) = 2;
                    else
                        dist = DistToThresh(Outer_upper_boundary,Outer_lower_boundary,ones_sum(pick_index,global_clock),current_time,Test_Limit);
                        costs(pick_index,(global_clock + 1)) = dist + (lambda_range(lambda_index) * current_time);
                    end
                end
            end
            
            if find(decisions == good_bin,1,'first')
                good_at = find(decisions == good_bin,1,'first');
                Trail_local_samples_to_find_1(trial_index) = sum(local_clocks(good_at,:));
                Trail_global_samples_to_find_1(trial_index) = sum(sum(
                Trial_samples_per_channel(trial_index) = mean(local_clocks(local_clocks ~= local_clocks(good_at)));
                
                if p_act_bin(good_at) ~= good_bin
                    Trial_error(trial_index) = 1;
                end
                break;
            end
            
            dummy = zeros(1,Histogram_max);
            for h = 1:1:Histogram_max
                dummy(h) = sum(local_clocks == h);
            end
            Trial_histogram(trial_index,:) = dummy;
        end
        %% Result Collection
        Mean_local_samples_to_find_1(channel_index,lambda_index) = mean(Trail_local_samples_to_find_1);
        Mean_global_samples_to_find_1(channel_index,lambda_index) = mean(Trail_global_samples_to_find_1);
        Mean_error(channel_index,lambda_index) = mean(Trial_error);
        Mean_samples_per_channel(channel_index,lambda_index) = mean(Trial_samples_per_channel);
        for h = 1:1:Histogram_max
            histogram(channel_index,lambda_index,h) = mean(Trial_histogram(:,h));
        end
    end
end


matlabpool close;

%% Save data in file
%store the results for comparison
fname = strcat('VarLambdaData-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);