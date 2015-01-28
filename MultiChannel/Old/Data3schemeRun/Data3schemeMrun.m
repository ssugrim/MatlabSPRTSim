%Generates "Data" for a 3 scheme analysis

%% Paths and local functions

addpath '/home/ssugrim/My Documents/Academics/research/ChannelSurfing/Simulation/SampleGen/lib'

matlabpool open;

%for use during cmd line run in linux

%Print diagonstic Messages
DEBUG = false;

%% GLOBAL PARAMETERS  

%Our sample allotment
deadline = 1000;

%number of Trials to test
Trials = 10000;

%Histogram bins
Histo_bin = 0:1:30;

%Number of Channels to test aganst
Channels = 2.^(5:1:10);
%Channels = [32];
%number of good channels we want to dig out (a preformance threhold)
good = 4;

%penalitly coefficent;
lambda = 0.05;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

delta = 0.1;
p_prime = 0.2;
split = 0.3;

%% DERVIED VALUES 

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%%%%%%%%%% right side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_rs = 1 - p_prime;

%boundaries of accept and reject region
p_1_rs = p_prime_rs + (delta * (1 - split));
p_0_rs = p_prime_rs - (delta * split);

%generate test lines
[L_1_rs,L_0_rs] = GenTestLines(alpha,beta,p_1_rs,p_0_rs,deadline);

%%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_ls = p_prime;

%boundaries of accept and reject region
p_1_ls = p_prime_ls + (delta * (1 - split));
p_0_ls = p_prime_ls - (delta * split);

%generate test lines
[L_1_ls,L_0_ls] = GenTestLines(alpha,beta,p_1_ls,p_0_ls,deadline);

%%% Convienence Aliases %%%

%Upper threshold
T_U = L_1_rs;

%Lower threshold
T_L = L_0_ls;

%Center Upper threshold
T_C_U = L_0_rs;

%Center Lower threshold
T_C_L = L_1_ls;

%maximum cost
max_cost =  pdist([deadline,T_U(deadline);deadline,T_L(deadline)],'euclidean') + (lambda * deadline); 

%% DATA STORAGE VARS - We do analysis on these

%Time to find the required number of good channels
Data_good_decision_time = zeros(length(Channels),Trials);

%Number of good channels available
Data_good_available = zeros(length(Channels),Trials);

%% SPRT Data

%Number of channels that were Significant to SPRT  - completed to a decision
Data_SPRT_Significant = zeros(length(Channels),Trials);

%Errors made during the trial
Data_SPRT_error_1 = zeros(length(Channels),Trials);
Data_SPRT_error_2 = zeros(length(Channels),Trials);

%Mean number of measurements per channel for SPRT
Data_SPRT_mean_m_k = zeros(length(Channels),Trials);

% Percent Classified by SPRT
Data_SPRT_percent_classified = zeros(length(Channels),Trials);

%SPRT allocation histogram
Data_SPRT_histogram = zeros(length(Channels),Trials,length(Histo_bin));

%% Simple Data

%Number of channels that were Significant to Simple  - At least one measurement
Data_Simple_Significant = zeros(length(Channels),Trials);

%Errors made with the simple Stratgey
Data_Simple_error_1 = zeros(length(Channels),Trials);
Data_Simple_error_2 = zeros(length(Channels),Trials);

%Mean number of measurements per channel for Simple
Data_Simple_mean_m_k = zeros(length(Channels),Trials);

% Percent Classified by Simple
Data_Simple_percent_classified = zeros(length(Channels),Trials);

%Simple allocation histogram
Data_Simple_histogram = zeros(length(Channels),Trials,length(Histo_bin));

%% Tree Data

%Number of channels that were Significant to Tree  - More than one pass
Data_Tree_Significant = zeros(length(Channels),Trials);

%Errors made with the simple Stratgey
Data_Tree_error_1 = zeros(length(Channels),Trials);
Data_Tree_error_2 = zeros(length(Channels),Trials);

%Mean number of measurements per channel for TREE
Data_Tree_mean_m_k = zeros(length(Channels),Trials);

% Percent Classified by Simple
Data_Tree_percent_classified = zeros(length(Channels),Trials);

%Simple allocation histogram
Data_Tree_histogram = zeros(length(Channels),Trials,length(Histo_bin));

%% Simulation Begins here 

%%%% Begin Trial Loop (i is my trial index) %%%%
for chan_num_ind = 1:1:length(Channels)
    K = Channels(chan_num_ind);   
    sprintf('On channel set %d',K)
    parfor i = 1:Trials
        %% cheap progress indicator
        if (mod(i,25) == 0)
            sprintf('On trail %d',i)
        end
        
        %% TRIAL SETUP - Per trial computed vales
        
        %per trial channel Parameters (rand draws from uniform 0,1)
        p_act = rand(1,K);
        
        %place the P_actuals in bins (the correct ones)
        bin_p_act = zeros(1,K);
        
        %bin the actual
        for k = 1:1:K
            if p_act(k) > p_prime_rs
                bin_p_act(k) = 3;
            elseif p_act(k) < p_prime_ls
                bin_p_act(k) = 1;
            else
                bin_p_act(k) = 2;
            end
        end
        
        %count the number of good one we could find on this trial
        Data_good_available(chan_num_ind,i) = sum(p_act < p_prime_ls);
        
        %% SPRT INITILIZATION
        
        dbg_print('Start SPRT',DEBUG);
        
        %cost matrix for each channel at each time instant j
        Trial_cost = zeros(K,deadline);
        
        %path taken by each channel
        Trial_path = zeros(K,deadline);
        
        %distance history
        Trial_dist = zeros(K,deadline);
        
        %indicies where good channels were found
        Trial_good_decision = zeros(1,good);
        
        %The first minimum out threshold distance
        Trial_e_k = zeros(1,K);
        
        %the current cost
        Trial_cur_cost = zeros(1,K);
        
        %What decision was made for channel K
        Trial_SPRT_decision = zeros(1,K);
        
        %when was it made (what global index)
        Trial_decided_at = zeros(1,K);
        
        %The sequence of Channels that was picked
        Trial_pick_sequence = zeros(1,deadline);
        
        %running sum of ones
        Trial_sprt_d_k = zeros(1,K);
        
        %number of measurements each channel k has recieved
        Trial_sprt_m_k = ones(1,K);
        
        %initalize cost and threshold distance
        for k = 1:1:K
            Trial_e_k(k) = min(abs(T_U(1) - Trial_sprt_d_k(k)),abs(T_L(1) - Trial_sprt_d_k(k)));
            Trial_cur_cost(k) = Trial_e_k(k);
        end
        
        %% SPRT Simulation LOOP - The actual SPRT run
        for j = 1:1:deadline+1
            %if All decision have been made stop this run
            if sum(Trial_SPRT_decision > 0 ) == K
                break;
            end
            
            %update the cost
            for k = 1:1:K
                if Trial_SPRT_decision(k) ~= 0
                    %if a decision was made cost is infinity, this channel needs
                    %no more measurements
                    Trial_cur_cost(k) = max_cost;
                else
                    %update the distance and cost
                    U_dist = pdist([Trial_sprt_m_k(k), T_U(Trial_sprt_m_k(k)); Trial_sprt_m_k(k),Trial_sprt_d_k(k)],'euclidean');
                    L_dist = pdist([Trial_sprt_m_k(k),T_L(Trial_sprt_m_k(k)); Trial_sprt_m_k(k),Trial_sprt_d_k(k)],'euclidean');
                    Trial_e_k(k) = min(U_dist,L_dist);
                    Trial_cur_cost(k) = Trial_e_k(k) + (lambda * Trial_sprt_m_k(k)) ;
                end
            end
            
            for k = 1:1:K
                %check for decisions
                if Trial_SPRT_decision(k) == 0
                    if Trial_sprt_d_k(k) > T_U(Trial_sprt_m_k(k))
                        %above the upper boundary is region 3
                        Trial_SPRT_decision(k) = 3;
                        Trial_decided_at(k) = j;
                    elseif Trial_sprt_d_k(k) < T_L(Trial_sprt_m_k(k))
                        %below the lower boundary is region 1, the "good"
                        %region since this means p is low
                        Trial_SPRT_decision(k) = 1;
                        Trial_decided_at(k) = j;
                        Trial_good_decision(find(Trial_good_decision == 0, 1, 'first')) = j;
                    elseif T_C_L(Trial_sprt_m_k(k)) < Trial_sprt_d_k(k) && Trial_sprt_d_k(k) < T_C_U(Trial_sprt_m_k(k))
                        %In between the center boundaries means we're in the
                        %second region
                        Trial_SPRT_decision(k) = 2;
                        Trial_decided_at(k) = j;
                    end
                end
            end
            
            %pick the channel smallest distance
            cur_pick = pick(Trial_cur_cost);
            Trial_pick_sequence(j) = cur_pick;
            %sample that channel
            Trial_sprt_d_k(cur_pick) = Trial_sprt_d_k(cur_pick) + random('bino',1,p_act(cur_pick));
            %increment it's meausrement count
            Trial_sprt_m_k(cur_pick) = Trial_sprt_m_k(cur_pick) + 1;
            
            %Store the trials cost vector (across all channels) on per indexed by time
            Trial_cost(:,j) = Trial_cur_cost';
            %Store the trials path vector (across all channels) on per indexed by time
            Trial_path(:,j) = Trial_sprt_d_k';
            %Store the trials minimum distance vector (across all channels) on per indexed by time
            Trial_dist(:,j) = Trial_e_k';
        end
        
        %% SPRT data Analysis - Collecting error and decision information
        
        % Use the bins to compute the error rates
        [Data_SPRT_error_1(chan_num_ind,i),Data_SPRT_error_2(chan_num_ind,i)] = TrialAlphaBeta(K,bin_p_act,Trial_SPRT_decision);
        dbg_str = sprintf('Trial SPRT error rates are alpha = %f, beta = %f ',Data_SPRT_error_1(chan_num_ind,i),Data_SPRT_error_2(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        %number of channels completed for trial i
        Data_SPRT_Significant(chan_num_ind,i) =  sum(Trial_SPRT_decision > 0);
        dbg_str = sprintf('SPRT charaterized %d at deadline, they were:',Data_SPRT_Significant(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        % Percent Classified by SPRT
        Data_SPRT_percent_classified(chan_num_ind,i) = Data_SPRT_Significant(chan_num_ind,i) / K;
        
        Trial_SPRT_p_hat = zeros(1,K); 
        %Compute mean p_hat
        for k = 1:1:K
            if Trial_sprt_m_k(k) ~= 0
               Trial_SPRT_p_hat(k) =  Trial_sprt_d_k(k) / (k);
            end
        end
        
        Data_SPRT_mean_p_hat_good(chan_num_ind,i) = mean(Trial_SPRT_p_hat(Trial_SPRT_decision == 1));
        
        %Store and report Mean number of samples used
        
        %Mean across all channels (including the zeros) Cuz they don't
        %understand what mean across measured means
        Data_SPRT_mean_m_k(chan_num_ind,i) = mean(Trial_sprt_m_k);
        dbg_str = sprintf('Mean number of samples used = %f ',Data_SPRT_mean_m_k(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        %how long did it take to find 4 good ones.
        Data_good_decision_time(chan_num_ind,i) = max(Trial_good_decision);
        
        %compute the SPRT histogram
        Data_SPRT_histogram(chan_num_ind,i,:) = TrialHistogram(Trial_sprt_m_k,Histo_bin);
        
        dbg_print('End SPRT',DEBUG);
                      
        %% Simple Scheme Initilization
        dbg_print('Start Simple',DEBUG);
        
        %running sum of ones
        Trial_simple_d_k = zeros(1,K);
        
        %number of measurements each channel k has recieved
        Trial_simple_m_k = zeros(1,K);
        
        %Simple scheme estimate of p_k
        Trial_simple_hat_p_k = zeros(1,K);
        
        %bins for the simple scheme
        Trial_simple_bin_hat_p = zeros(1,K);
        
        %% Simple Scheme Simulation Loop
        
        %scan from lowest to highest (taking count of how many measuments each
        %channel gets.
        for l = 1:1:deadline+1
            %Index wrap around
            cur_index = mod(l,K) + 1;
            
            %take a single sample
            Trial_simple_d_k(cur_index) = Trial_simple_d_k(cur_index) + random('bino',1,p_act(cur_index));
            
            %Count how many times this channel has been measured
            Trial_simple_m_k(cur_index) = Trial_simple_m_k(cur_index) + 1;
        end
        
        %% Simple Scheme Data Analysis
        
        %computes the estimates we have measrements for
        for k = 1:1:K
            if Trial_simple_m_k(k) ~= 0
                Trial_simple_hat_p_k(k) = Trial_simple_d_k(k) / Trial_simple_m_k(k);
            end
        end
        
        %bin the Estimates
        for k = 1:1:K
            if Trial_simple_hat_p_k(k) > p_prime_rs
                Trial_simple_bin_hat_p(k) = 3;
            elseif Trial_simple_hat_p_k(k) < p_prime_ls
                Trial_simple_bin_hat_p(k) = 1;
            else
                %if it was never measured, throw it in 2.
                Trial_simple_bin_hat_p(k) = 2;
            end
        end
        
        %Use the bins to compute the error rates
        [Data_Simple_error_1(chan_num_ind,i),Data_Simple_error_2(chan_num_ind,i)] = TrialAlphaBeta(K,bin_p_act,Trial_simple_bin_hat_p);
        dbg_str = sprintf('Trial Simple error rates are alpha = %f, beta = %f ',Data_Simple_error_1(chan_num_ind,i),Data_Simple_error_2(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        %Number of signifcant channels
        Data_Simple_Significant(chan_num_ind,i) = sum(Trial_simple_m_k > 0);
        dbg_str = sprintf('Number of measured channels = %d', Data_Simple_Significant(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        % Percent Classified by SPRT
        Data_Simple_percent_classified(chan_num_ind,i) = Data_Simple_Significant(chan_num_ind,i) / K;
        
        %Compute mean p_hat
        Data_Simple_mean_p_hat_good(chan_num_ind,i) = mean(Trial_simple_hat_p_k(Trial_simple_bin_hat_p == 1));
                
        %Store and report Mean number of samples used
        Data_Simple_mean_m_k(chan_num_ind,i) = mean(Trial_simple_m_k);
        dbg_str = sprintf('Mean number of samples used = %f ',Data_Simple_mean_m_k(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        %compute the Simple histogram
        Data_Simple_histogram(chan_num_ind,i,:) = TrialHistogram(Trial_simple_m_k,Histo_bin);
                
        dbg_print('End Simple',DEBUG);
        
        %% Tree Scheme Initilization
        dbg_print('Start Tree',DEBUG);
        
        %fraction of samples to keep in reserve for making decisons
        reserve = 0.5;
        
        %running sum of ones
        Trial_tree_d_k = zeros(1,K);
        
        %number of measurements each channel k has recieved
        Trial_tree_m_k = zeros(1,K);
        
        %slack factor - number of ones we are willing to miss
        slack = 0.25;
        
        %Samples per channel per pass
        pass_sample = 4;
        
        %% Tree Scheme Simulation Loop
        
        %scan from lowest to highest (taking count of how many measuments each
        %channel gets). Here we stop at deadline/2 since we're doing two
        %measurements per channel
        
        

        %do a first pass with using up the ((1 - reserve) * limit)
        %allotment
        
        used =  1;
        allotment = floor(((1 - reserve) * deadline) / pass_sample);
        
        for k = 1:1:min(allotment,K)
            for n = 1:1:pass_sample
                Trial_tree_d_k(k) = Trial_tree_d_k(k) + random('bino',1,p_act(k));
                Trial_tree_m_k(k) = Trial_tree_m_k(k) + 1;
                used = used + 1;
            end
        end
        
        %sprintf('Completed first pass, used = %d',used)
        
        %and now distribute measurements that are on the edges
        cur_index = 1;
        pass = 0;
        while used < deadline + 1
            %relaxed too much I must be done (have to let it get past .5
            %otherwise we miss some channels even though we have
            %measurements
            if slack > 0.6
                break;
            end
            
            picked = zeros(1,K);
            %pick a set of channels
            for k = 1:1:K
                if Trial_tree_m_k(k) > 0
                    ratio = Trial_tree_d_k(k) / Trial_tree_m_k(k);
                    % parameters are within some bound and the channels haven't
                    % passed some sane limit
                    if (ratio < slack || ratio > 1 - slack) && Trial_tree_m_k(k) < 50
                        picked(k) = 1;
                    end
                end
            end
            
            %Didn't find any thing worth picking, relax the criteria a bit.
            if sum(picked) == 0
                slack = slack + 0.05;
                %sprintf('Bumped up the slack to %f', slack)
                continue;
            end
            
            %sprintf('Pass = %d, # picked = %d',pass, length(picked(picked > 0)))
            
            %distribute measurements to the picks
            for k = 1:1:K
                if picked(k) ~= 0
                    for n = 1:1:pass_sample
                        Trial_tree_d_k(k) = Trial_tree_d_k(k) + random('bino',1,p_act(k));
                        Trial_tree_m_k(k) = Trial_tree_m_k(k) + 1;
                        used = used + 1;
                    end
                end
            end
            
            %didn't use any measurements so I must be done
            pass = pass + 1;
        end
        
        %% Tree Scheme Data Analysis
        
        %Simple scheme estimate of p_k
        Trial_tree_hat_p_k = zeros(1,K);
        
        %computes the estimates of paramters for channels we have measrements for
        for k = 1:1:K
            if Trial_tree_m_k(k) > pass_sample
                Trial_tree_hat_p_k(k) = Trial_tree_d_k(k) / Trial_tree_m_k(k);
            else
                Trial_tree_hat_p_k(k) = -1;
            end
        end
        
        %bins for the simple scheme
        Trial_tree_bin_hat_p = zeros(1,K);
        
        %bin the Estimates
        for k = 1:1:K
            if Trial_tree_hat_p_k(k) ~= -1
                if Trial_tree_hat_p_k(k) > p_prime_rs
                    Trial_tree_bin_hat_p(k) = 3;
                elseif Trial_tree_hat_p_k(k) < p_prime_ls
                    Trial_tree_bin_hat_p(k) = 1;
                else
                    %if it was never measured, throw it in 2.
                    Trial_tree_bin_hat_p(k) = 2;
                end
            end
        end
        
        %Use the bins to compute the error rate
        [Data_Tree_error_1(chan_num_ind,i),Data_Tree_error_2(chan_num_ind,i)] = TrialAlphaBeta(K,bin_p_act, Trial_tree_bin_hat_p);
        dbg_str = sprintf('Trial Tree error rates are alpha = %f, beta = %f ',Data_Tree_error_1(chan_num_ind,i),Data_Tree_error_2(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        %number of channels with a significant measurement
        Data_Tree_Significant(chan_num_ind,i) = length(Trial_tree_m_k(Trial_tree_m_k > pass_sample));
        dbg_str = sprintf('Number of channels with significant= %d', Data_Tree_Significant(chan_num_ind,i) );
        dbg_print(dbg_str,DEBUG);
        
        % Percent Classified by SPRT
        Data_Tree_percent_classified(chan_num_ind,i) = Data_Tree_Significant(chan_num_ind,i) / K;
           
        %Store and report Mean number of samples used
        Data_Tree_mean_m_k(chan_num_ind,i) = mean(Trial_tree_m_k);
        dbg_str = sprintf('Mean number of samples used = %f ',Data_Tree_mean_m_k(chan_num_ind,i));
        dbg_print(dbg_str,DEBUG);
        
        %compute the Tree histogram
        Data_Tree_histogram(chan_num_ind,i,:) = TrialHistogram(Trial_tree_m_k,Histo_bin);
       
        dbg_print('End Tree',DEBUG);
        
        %% Sanity Diagonstic output
        if DEBUG
            %         sprintf('Bins')
           % bin_p_act
           % Trial_SPRT_decision
            %Trial_simple_bin_hat_p
            %Trial_tree_bin_hat_p
            %
            %         sprintf('Estimates')
           % p_act
            %         Trial_sprt_m_k
            %         Trial_simple_hat_p_k
            %        Trial_simple_d_k
            %        Trial_simple_m_k
            %        Trial_tree_hat_p_k
            %        Trial_tree_d_k
            %        Trial_tree_m_k
        end
        
        %% End Trial
    end
end
%% Store the results for comparison
matlabpool close;
fname = sprintf('Data3SchemeTri%dmin%dmax%d.mat',Trials,min(Channels),max(Channels));
save(fname);
exit;
