%Computes the number of samples required to find 1 good channel in a set of
%varying size with varying lambda

%for use during cmd line run in linux
%addpath(genpath('/home/ssugrim/My Documents/Academics/research/ChannelSurfing/Simulation/SampleGen/'))

%Multi Threading
matlabpool open;

%% Test Tuning Params
%global Parameters

debug = false; 

%number of channels
Channels = [2^10];

%The good channel value
good_val = 0.03;

%The bin we deem as the desriable one
good_bin = 1;

%penalitly coefficent;
lambda_range = 0:.005:0.95;

%The hard limt All tests must stop at (Lenght of our vectors, should be
%longer than eta, for some allowance).
Test_Limit = 1000;

Histogram_max = 100;

%number of Trials
Trial_limit = 500;

%test line recentering  fudge factor
fudge_factor = 0.95;

%forgiveness fudge factor to prevent early failures
%how many mistakes we're willing to tolerate
forgive_allow = 2;

%Clock threshold after which we don't forgive. 
forgive_clock = 10;

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
Outer_lower_boundary = L_0_ls - fudge_factor;
Outer_upper_boundary = L_1_rs - fudge_factor;
Inner_lower_boundary = L_1_ls - fudge_factor;
Inner_upper_boundary = L_0_rs - fudge_factor;

%% Storage Vars for plot generation

Mean_local_samples_to_find_1 = zeros(length(Channels),length(lambda_range));
Mean_global_samples_to_find_1 = zeros(length(Channels),length(lambda_range));
Mean_error = zeros(length(Channels),length(lambda_range));
Mean_channels_searched = zeros(length(Channels),length(lambda_range));
Mean_samples_per_channel = zeros(length(Channels),length(lambda_range));
Mean_forgivness_given = zeros(length(Channels),length(lambda_range));
Mean_value_found_estimate = zeros(length(Channels),length(lambda_range));
Mean_value_found_actual = zeros(length(Channels),length(lambda_range));
Mean_trial_failure = zeros(length(Channels),length(lambda_range));

% if debug
%     %% Global debug storage

%     Mean_estimate_of_good = zeros(length(Channels),length(lambda_range));
%     histogram = zeros(length(Channels),length(lambda_range),Histogram_max);
% end

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
        
        %number of chanels searched
        Trial_channels_searched = zeros(1,Trial_limit);
        
        %Mean number of sammples perchannel during this trial
        Trial_samples_per_channel = zeros(1,Trial_limit);
        
        %number of times forgivness was given
        Trial_forgivness_given = zeros(1,Trial_limit);
        
        % Estimate of the value found
        Trial_value_found_estimate = zeros(1,Trial_limit);
        
        %What the value actually was
        Trial_value_found_actual = zeros(1,Trial_limit);
        
        %was the trial a failure
        Trial_failure = zeros(1,Trial_limit);
            
        
%         if debug
%             %% Trial Debug Storage
%             %parameter estimate for the channel we found
%             Trial_estimate_of_good = zeros(1,Trial_limit);
%             
%             %parameter estimate for the channel we found
%             Trial_index_path = zeros(Trial_limit,Test_Limit);
%             
%             %actualy values of P
%             Trial_p_act = zeros(Trial_limit,Channels(channel_index));
%             
%             %The local clock for each channel per trail
%             Trial_local_clocks = zeros(Trial_limit,Channels(channel_index));
%             
%             %The ones sum for each channel per trail
%             Trial_ones_sum = zeros(Trial_limit,Channels(channel_index));
%             
%             %Decsions made if any
%             Trial_decisions = zeros(Trial_limit,Channels(channel_index));
%             
%             %parameter estimate for the channel we found
%             Trial_good_act = zeros(1,Trial_limit);
%             
%             %The histogram for this trial
%             Trial_histogram = zeros(Trial_limit,Histogram_max);
%             
%             %Samples for plotting
%             Trial_sample_path = zeros(Trial_limit,Channels(channel_index),Test_Limit);
%             
%             %Cost for debugging
%             Trial_cost_path = zeros(Trial_limit,Channels(channel_index),Test_Limit);
%         end
        
        %% computed Trail Values
        %maximal Right now, just some really large number so I don't come
        %back
        max_cost = 50000;
        
        parfor trial_index = 1:Trial_limit
        %for trial_index = 1:1:Trial_limit
            %% Trial Loop
            sprintf('On Trial %d',trial_index)
            
            %% Genrerate a Random vector of appropriate length then bin it
            p_act = rand(1,Channels(channel_index));
            
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
            
            %% Trail temprorary vars
            
            %One clock per channel - Number of measurements each channel
            %gets (m_k)
            local_clocks = zeros(1,Channels(channel_index));
            
            %Vector of costs for each channel (e_k)
            costs = zeros(1,Channels(channel_index));
            
            %vector of decision for each channel, zero if undecided, 1,2 or
            %3 for the respective bins
            decisions = zeros(1,Channels(channel_index));
            
            %vector of accumulated ones - For this set of parameters the
            %fair starting point is (1,1).
            ones_sum = zeros(1,Channels(channel_index));
            
            %% Trial Computed Values
            
            %Compute the initial cost for each channel (should all be the
            %same). Since we're centering on 1,1 there should be no
            %"measurement penality"
            for k = 1:1:Channels(channel_index)
               % costs(k) = min(pdist([1,ones_sum(1);1,Outer_lower_boundary(1)]),pdist([1,ones_sum(1);1,Outer_upper_boundary(1)]));
                costs(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,ones_sum(k),local_clocks(k),Test_Limit);
            end
            
            %path storage Vars
%             index_path = zeros(1,Test_Limit);
%             sample_path = zeros(Channels(channel_index),Test_Limit);
%             cost_path = zeros(Channels(channel_index),Test_Limit);
            
            for global_clock = 1:1:Test_Limit
                %% The actual SPRT Loop
                
                % Check if the test has finished - It's finished when we
                % pick a channel and place it in the good bin. 
                if find(decisions == good_bin,1,'first')
                   break;
                end
                
 %               cost_path(:,global_clock) = costs;
                
                %pick a channel
                pick_index =  pick(costs);
                                
  %              index_path(global_clock) = pick_index;
                
                %update that channels local clock
                local_clocks(pick_index) = local_clocks(pick_index) + 1;
                %sprintf('Update on chan:%d, p_act:%f, cost:%f, ones:%d, LeftT:%f, RightT%f, lclock:%d',pick_index, p_act(pick_index), costs(pick_index),ones_sum(pick_index),Outer_upper_boundary(local_clocks(pick_index)),Outer_lower_boundary(local_clocks(pick_index)),local_clocks(pick_index))
                              
                %Sample that channel and add the result to the running ones
                ones_sum(pick_index) = ones_sum(pick_index) +  random('bino',1,p_act(pick_index));
                %sprintf('Update on chan:%d, p_act:%f, cost:%f, ones:%d, LeftT:%f, RightT%f, lclock:%d',pick_index, p_act(pick_index), costs(pick_index),ones_sum(pick_index),Outer_upper_boundary(local_clocks(pick_index)),Outer_lower_boundary(local_clocks(pick_index)),local_clocks(pick_index))
%                sample_path(pick_index,local_clocks(pick_index)) = ones_sum(pick_index);
               
                %update the cost and %check for decisons 
                for k = 1:1:Channels(channel_index)
                    if  decisions(k) == 0
                        if ones_sum(k) > Outer_upper_boundary(local_clocks(k)+1)
                            %decide region 3
                %            sprintf('Dec:3,chan:%d,sum:%d,bound:%f,act:%f,mea:%d,cost:%f',k,ones_sum(k),Outer_upper_boundary(local_clocks(k)),p_act(k),local_clocks(k),costs(k))
                            costs(k) = max_cost;
                            decisions(k) = 3;
                        elseif ones_sum(k) < Outer_lower_boundary(local_clocks(k)+1)
                            %decide region 1
                      %      sprintf('Dec:1,chan:%d,sum:%d,bound:%f,act:%f,mea:%d,cost:%f',k,ones_sum(k),Outer_lower_boundary(local_clocks(k)),p_act(k),local_clocks(k),costs(k))
                            costs(k) = max_cost;
                            decisions(k) = 1;
                        elseif Inner_lower_boundary(local_clocks(k)+1) < ones_sum(k) && ones_sum(k) < Inner_upper_boundary(local_clocks(k)+1)
                            %decide region 2
                       %     sprintf('Dec:2,chan:%d,sum:%d,bound:%f_%f,act:%f,mea:%d,cost:%f',k,ones_sum(k),Inner_lower_boundary(local_clocks(k)),Inner_upper_boundary(local_clocks(k)),p_act(k),local_clocks(k),costs(k))
                            costs(k) = max_cost;
                            decisions(k) = 2;
                        else
                            %keep sampling
                          %  u_dist = pdist([local_clocks(k),ones_sum(k);local_clocks(k),Outer_lower_boundary(local_clocks(k))]);
                          %  l_dist = pdist([local_clocks(k),ones_sum(k);local_clocks(k),Outer_upper_boundary(local_clocks(k))]);
                          %  costs(k) = min(u_dist,l_dist) + (lambda_range(lambda_index) * local_clocks(k));
                          if forgive_clock > local_clocks(k)
                              % if we are have only made 2 mistakes above
                              % or below
                              if (forgive_allow > ones_sum(k)) || ((local_clocks(k) - ones_sum(k)) < forgive_allow)
                              forgivness = ForgiveAmt(Outer_upper_boundary,Outer_lower_boundary,ones_sum(k),local_clocks(k),Test_Limit)
                              %Track the ammount of forgiveness used during a trail.
                              Trial_forgivness_given(trial_index) = Trial_forgivness_given(trial_index) + 1;
                              end
                          else
                              forgivness = 0;
                          end
                            costs(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,ones_sum(k),local_clocks(k),Test_Limit) + (local_clocks(k) * lambda_range(lambda_index)) + forgivness;
                        end
                    end
                end
            end
            
            %Incase we didn't finish
            good_at = find(decisions == good_bin,1,'first');
            if good_at
                %How many samples did the good channel take?
                Trail_local_samples_to_find_1(trial_index) = local_clocks(good_at);
                %was a error made?
                if p_act_bin(good_at) ~= good_bin
                    Trial_error(trial_index) = 1;
                end
                
                %How man channels did we have to go through before this
                %one? Any thing with a non zero clock - 1 cuz we found one
                bad_clocks = local_clocks([1:(good_at - 1), (good_at + 1):end]);
                Trial_channels_searched(trial_index) = sum(bad_clocks(bad_clocks ~= 0));
                
                %Mean number of sammples perchannel during this trial
                Trial_samples_per_channel(trial_index) = mean(bad_clocks(bad_clocks ~= 0));
                
                
                % Estimate of the value found
                Trial_value_found_estimate(trial_index) = ones_sum(good_at) / local_clocks(good_at);
                
                %What the value actually was
                Trial_value_found_actual(trial_index) = p_act(good_at);
                
            else
                %Nonsense number
                Trial_failure(trial_index) = 1;
                Trail_local_samples_to_find_1(trial_index) = Test_Limit;
                Trial_channels_searched(trial_index) = sum(local_clocks(local_clocks ~= 0));
               % Trial_samples_per_channel(trial_index) = mean(local_clocks(local_clocks ~= 0));
            end
            
            sprintf('Finished trial, lambda %f, local:%d, global:%d',lambda_range(lambda_index), local_clocks(good_at),global_clock)
            
            Trail_global_samples_to_find_1(trial_index) = global_clock;
            
%             if debug
%                 %% for debug purposes
%                 Trial_samples_per_channel(trial_index) = mean(local_clocks(local_clocks ~= local_clocks(good_at)));
%                 Trial_channels_searched(trial_index) = sum(local_clocks ~= 0);
%                 Trial_estimate_of_good(trial_index) = ones_sum(good_at) / local_clocks(good_at);
%                 Trial_index_path(trial_index,:) = index_path;
%                 Trial_p_act(trial_index,:) = p_act;
%                 Trial_local_clocks(trial_index,:) = local_clocks;
%                 Trial_ones_sum(trial_index,:) = ones_sum;
%                 Trial_decisions(trial_index,:) = decisions;
%                 Trial_good_act(trial_index) = p_act(good_at);
%                 
%                 dummy = zeros(1,Histogram_max);
%                 for h = 1:1:Histogram_max
%                     dummy(h) = sum(local_clocks == h - 1);
%                 end
%                 Trial_histogram(trial_index,:) = dummy;
%                 
%                 for k = 1:Channels(channel_index)
%                     Trial_sample_path(trial_index,k,:) = sample_path(k,:);
%                     Trial_cost_path(trial_index,k,:) = cost_path(k,:);
%                 end
%                 
%             end
        end
        %% Result Collection
        %We subtract one for each of these because we started at (1,1) to
        %make the test fair, thus we overcounted by 1 sample. 
        Mean_local_samples_to_find_1(channel_index,lambda_index) = mean(Trail_local_samples_to_find_1);
        Mean_global_samples_to_find_1(channel_index,lambda_index) = mean(Trail_global_samples_to_find_1);
        Mean_error(channel_index,lambda_index) = mean(Trial_error);
        Mean_channels_searched(channel_index,lambda_index) = mean(Trial_channels_searched);
        Mean_samples_per_channel(channel_index,lambda_index) = mean(Trial_samples_per_channel);
        Mean_forgivness_given(channel_index,lambda_index) = mean(Trial_forgivness_given);
        Mean_value_found_estimate(channel_index,lambda_index) = mean(Trial_value_found_estimate);
        Mean_value_found_actual(channel_index,lambda_index) = mean(Trial_value_found_actual);
        Mean_trial_failure(channel_index,lambda_index) = mean(Trial_failure);
        
        
%         %% Debug Globals
%         Mean_channels_searched(channel_index,lambda_index) = mean(Trial_channels_searched);
%         Mean_estimate_of_good(channel_index,lambda_index) = mean(Trial_estimate_of_good);
%         Mean_samples_per_channel(channel_index,lambda_index) = mean(Trial_samples_per_channel);
%         for h = 1:1:Histogram_max
%             histogram(channel_index,lambda_index,h) = mean(Trial_histogram(:,h));
%         end
    end
end


matlabpool close;

%% Save data in file
%store the results for comparison
fname = strcat('VarLambdaDataNewDist-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);
