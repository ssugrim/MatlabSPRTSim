%Generates "Data" for a Multiple parameter choices of the sprt

%% Paths and local functions

matlabpool open;

%for use during cmd line run in linux
%addpath '../lib'

%Print diagonstic Messages
DEBUG = false;

%% GLOBAL PARAMETERS  

%Our sample allotment
deadline = 1000;

%number of Trials to test
Trials = 5;

%Histogram bins
Histo_bin = 0:1:30;

%Number of Channels to test aganst
Channels = 2.^(5:1:10);

%% SPRT parameter sets

%penalitly coefficent;
lambda = 0.05;

%alpha and beta design choices
alpha = 0.01;
beta = 0.05;

delta_range = (0.15);
p_prime = 0.2;
split = 0.5;

%% DATA STORAGE VARS - We do analysis on these

%Number of good channels available
Data_good_available = zeros(length(delta_range),Trials,length(Channels));

%% SPRT Data

%Number of channels that were Significant to SPRT  - completed to a decision
Data_SPRT_Significant = zeros(length(delta_range),Trials,length(Channels));

%Errors made during the trial
Data_SPRT_error_1 = zeros(length(delta_range),Trials,length(Channels));
Data_SPRT_error_2 = zeros(length(delta_range),Trials,length(Channels));

%Mean number of measurements per channel for SPRT
Data_SPRT_mean_m_k = zeros(length(delta_range),Trials,length(Channels));

% Percent Classified by SPRT
Data_SPRT_percent_classified = zeros(length(delta_range),Trials,length(Channels));

% Percent Average Bin Value
Data_SPRT_mean_bin_1 = zeros(length(delta_range),Trials,length(Channels));
Data_SPRT_mean_bin_2 = zeros(length(delta_range),Trials,length(Channels));
Data_SPRT_mean_bin_3 = zeros(length(delta_range),Trials,length(Channels));

%Mean samples between decsion
Data_SPRT_mean_samples_between_decsion = zeros(length(delta_range),Trials,length(Channels));

%means samples to make a decision
Data_SPRT_mean_samples_make_decision = zeros(length(delta_range),Trials,length(Channels));

%SPRT allocation histogram
Data_SPRT_histogram = zeros(length(delta_range),Trials,length(Channels),length(Histo_bin));

%% DERVIED VALUES 

for delta_ind  = 1:length(delta_range)
    delta = delta_range(delta_ind);
    sprintf('On Delta %f',delta)
    
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
    
    %maximum cost
    max_cost =  pdist([deadline,T_U(deadline);deadline,T_L(deadline)],'euclidean') + (lambda * deadline);
    
    %% Simulation Begins here
    
    %%%% Begin Trial Loop (i is my trial index) %%%%
    for chan_num_ind = 1:length(Channels)
        K = Channels(chan_num_ind);
        sprintf('On channel set %d',K)
        parfor trial_ind = 1:Trials
       %  for trial_ind = 1:Trials
            %% cheap progress indicator
            if (mod(trial_ind,25) == 0)
                sprintf('On trail %d',trial_ind)
            end
            
            %% TRIAL SETUP - Per trial computed vales
            
            %per trial channel Parameters (rand draws from uniform 0,1)
            p_act = rand(1,K);
            
            %% SPRT INITILIZATION
            
            dbg_print('Start SPRT',DEBUG);
            
            %cost matrix for each channel at each time instant j
            Trial_cost = zeros(K,deadline);
            
            %path taken by each channel
            Trial_path = zeros(K,deadline);
            
            %distance history
            Trial_dist = zeros(K,deadline);
            
            %The minimum outer threshold distance
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
            
            %mean number of samples between decision;
            Trail__samples_between_decsion = zeros(1,K);
            Trial_samples_since_last_decision = 0;
            
            
            
            
            %initalize cost and threshold distance
            for k = 1:1:K
                Trial_e_k(k) = min(abs(T_U(1) - Trial_sprt_d_k(k)),abs(T_L(1) - Trial_sprt_d_k(k)));
                Trial_cur_cost(k) = Trial_e_k(k);
            end
            
            %% SPRT Simulation LOOP - The actual SPRT run
            for simu_index = 1:1:deadline+1
                %if All decision have been made stop this run
                if sum(Trial_SPRT_decision > 0 ) == K
                    break;
                end
                
                %bump up the sample count
                Trial_samples_since_last_decision = Trial_samples_since_last_decision + 1;
                
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
                            Trail__samples_between_decsion(k) = Trial_samples_since_last_decision;
                            Trial_samples_since_last_decision = 0;
                            Trial_decided_at(k) = simu_index;
                        elseif Trial_sprt_d_k(k) < T_L(Trial_sprt_m_k(k))
                            %below the lower boundary is region 1, the "good"
                            %region since this means p is low
                            Trial_SPRT_decision(k) = 1;
                            Trail__samples_between_decsion(k) = Trial_samples_since_last_decision;
                            Trial_samples_since_last_decision = 0;
                            Trial_decided_at(k) = simu_index;
                        elseif T_C_L(Trial_sprt_m_k(k)) < Trial_sprt_d_k(k) && Trial_sprt_d_k(k) < T_C_U(Trial_sprt_m_k(k))
                            %In between the center boundaries means we're in the
                            %second region
                            Trial_SPRT_decision(k) = 2;
                            Trail__samples_between_decsion(k) = Trial_samples_since_last_decision;
                            Trial_samples_since_last_decision = 0;
                            Trial_decided_at(k) = simu_index;
                        end
                    end
                end
                
                %pick the channel smallest distance
                cur_pick = pick(Trial_cur_cost);
                Trial_pick_sequence(simu_index) = cur_pick;
                %sample that channel
                Trial_sprt_d_k(cur_pick) = Trial_sprt_d_k(cur_pick) + random('bino',1,p_act(cur_pick));
                %increment it's meausrement count
                Trial_sprt_m_k(cur_pick) = Trial_sprt_m_k(cur_pick) + 1;
                
                %Store the trials cost vector (across all channels) on per indexed by time
                Trial_cost(:,simu_index) = Trial_cur_cost';
                %Store the trials path vector (across all channels) on per indexed by time
                Trial_path(:,simu_index) = Trial_sprt_d_k';
                %Store the trials minimum distance vector (across all channels) on per indexed by time
                Trial_dist(:,simu_index) = Trial_e_k';
            end
            
            %% SPRT data Analysis - Collecting error and decision information
            
            
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
            Data_good_available(delta_ind,trial_ind,chan_num_ind) = sum(p_act < p_prime_ls);
            
            
            % Use the bins to compute the error rates
            [Data_SPRT_error_1(delta_ind,trial_ind,chan_num_ind),Data_SPRT_error_2(delta_ind,trial_ind,chan_num_ind)] = TrialAlphaBeta(K,bin_p_act,Trial_SPRT_decision);
            dbg_str = sprintf('Trial SPRT error rates are alpha = %f, beta = %f ',Data_SPRT_error_1(delta_ind,trial_ind,chan_num_ind),Data_SPRT_error_2(delta_ind,trial_ind,chan_num_ind));
            dbg_print(dbg_str,DEBUG);
            
            %number of channels completed for trial i
            Data_SPRT_Significant(delta_ind,trial_ind,chan_num_ind) =  sum(Trial_SPRT_decision > 0);
            dbg_str = sprintf('SPRT charaterized %d at deadline, they were:',Data_SPRT_Significant(delta_ind,trial_ind,chan_num_ind));
            dbg_print(dbg_str,DEBUG);
            
            % Percent Classified by SPRT
            Data_SPRT_percent_classified(delta_ind,trial_ind,chan_num_ind) = Data_SPRT_Significant(delta_ind,trial_ind,chan_num_ind) / K;
            
            Trial_SPRT_p_hat = zeros(1,K);
            %Compute mean p_hat
            for k = 1:1:K
                if Trial_sprt_m_k(k) ~= 0 && Trial_sprt_d_k(k) ~= 0
                    Trial_SPRT_p_hat(k) =  Trial_sprt_d_k(k) / Trial_sprt_m_k(k);
                end
            end
            
            % Percent Average Bin Value
            Data_SPRT_mean_bin_1(delta_ind,trial_ind,chan_num_ind) = mean(Trial_SPRT_p_hat(Trial_SPRT_decision == 1));
            Data_SPRT_mean_bin_2(delta_ind,trial_ind,chan_num_ind) = mean(Trial_SPRT_p_hat(Trial_SPRT_decision == 2));
            Data_SPRT_mean_bin_3(delta_ind,trial_ind,chan_num_ind) = mean(Trial_SPRT_p_hat(Trial_SPRT_decision == 3));
            
            %Mean samples between decsion
            Data_SPRT_mean_samples_between_decsion(delta_ind,trial_ind,chan_num_ind) = mean(Trail__samples_between_decsion(Trail__samples_between_decsion ~= 0));
            
            %means samples to make a decision
            Data_SPRT_mean_samples_make_decision(delta_ind,trial_ind,chan_num_ind) = mean(Trial_sprt_m_k(Trial_SPRT_decision ~= 0));
            
            %Mean across all channels (including the zeros) Cuz they don't
            %understand what mean across measured means
            Data_SPRT_mean_m_k(delta_ind,trial_ind,chan_num_ind) = mean(Trial_sprt_m_k);
            dbg_str = sprintf('Mean number of samples used = %f ',Data_SPRT_mean_m_k(delta_ind,trial_ind,chan_num_ind));
            dbg_print(dbg_str,DEBUG);
            
            %compute the SPRT histogram
            Data_SPRT_histogram(delta_ind,trial_ind,chan_num_ind,:) = TrialHistogram(Trial_sprt_m_k,Histo_bin);
            
            dbg_print('End SPRT',DEBUG);
            
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
end
%% Store the results for comparison
matlabpool close;
fname = sprintf('Data3SchemeTri%dmin%dmax%d.mat',Trials,min(Channels),max(Channels));
save(fname);

