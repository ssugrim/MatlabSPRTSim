%Generates "Data" for a Multiple lambda choices of the sprt

%% Paths and local functions

matlabpool open;

%for use during cmd line run in linux
%addpath '../lib'

%Print diagonstic Messages
DEBUG = true;

%% GLOBAL PARAMETERS  

%Our sample allotment
Test_limit = 1000;

%number of Trials to test
Trials = 100;

%Histogram bins
Histo_bin = 0:1:30;

%Number of Channels to test aganst
Channels = 2^10;

%% SPRT parameter sets

%thresholds
delta = 0.125;
p_prime = 0.2;
split = 0.5;

%alpha and beta design choices
alpha = 0.01;
beta = 0.05;

%penalitly coefficent;
lambda_range = 0:0.05:0.95;

%test line recentering  fudge factor
fudge_factor = 0.95;

%forgiveness fudge factor to prevent early failures
%how many mistakes we're willing to tolerate
forgive_allow = 2;

%Clock threshold after which we don't forgive. 
forgive_clock = 10;


%% DATA STORAGE VARS - We do analysis on these


%% SPRT Data

%Number of channels that were Significant to SPRT  - completed to a decision
Data_SPRT_Significant = zeros(length(lambda_range),Trials,length(Channels));

%Errors made during the trial
Data_SPRT_error_1 = zeros(length(lambda_range),Trials,length(Channels));
Data_SPRT_error_2 = zeros(length(lambda_range),Trials,length(Channels));

%Mean number of measurements per channel for SPRT
Data_SPRT_mean_m_k = zeros(length(lambda_range),Trials,length(Channels));

% Fraction Classified by SPRT
Data_SPRT_fraction_classified = zeros(length(lambda_range),Trials,length(Channels));

% Percent Average Bin Value
Data_SPRT_mean_bin_1 = zeros(length(lambda_range),Trials,length(Channels));
Data_SPRT_mean_bin_2 = zeros(length(lambda_range),Trials,length(Channels));
Data_SPRT_mean_bin_3 = zeros(length(lambda_range),Trials,length(Channels));

%Mean samples between decsion
Data_SPRT_mean_samples_between_decsion = zeros(length(lambda_range),Trials,length(Channels));

%means samples to make a decision
Data_SPRT_mean_samples_make_decision = zeros(length(lambda_range),Trials,length(Channels));

%means samples to make a decision
Data_SPRT_mean_channels_passed_on = zeros(length(lambda_range),Trials,length(Channels));

%SPRT allocation histogram
Data_SPRT_histogram = zeros(length(lambda_range),Trials,length(Channels),length(Histo_bin));



for lambda_ind  = 1:length(lambda_range)
    lambda = lambda_range(lambda_ind);
    sprintf('On lambda %f',lambda)
    
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
    [L_1_rs,L_0_rs] = GenTestLines(alpha,beta,p_1_rs,p_0_rs,Test_limit);
    
    %%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%
    
    %tuneing parameters - center and width
    p_prime_ls = p_prime;
    
    %boundaries of accept and reject region
    p_1_ls = p_prime_ls + (delta * (1 - split));
    p_0_ls = p_prime_ls - (delta * split);
    
    %generate test lines
    [L_1_ls,L_0_ls] = GenTestLines(alpha,beta,p_1_ls,p_0_ls,Test_limit);
    
    %% Convienence Aliases
    %plus slight recentering
    Outer_lower_boundary = L_0_ls - fudge_factor;
    Outer_upper_boundary = L_1_rs - fudge_factor;
    Inner_lower_boundary = L_1_ls - fudge_factor;
    Inner_upper_boundary = L_0_rs - fudge_factor;
        
    %maximum cost - Just some large number
    max_cost =  50000;
    
    %% Simulation Begins here
    
    %%%% Begin Trial Loop (i is my trial index) %%%%
    for chan_num_ind = 1:length(Channels)
        N = Channels(chan_num_ind);
        
        sprintf('On channel set %d',N)
        parfor trial_ind = 1:Trials
       %for trial_ind = 1:Trials
            %% cheap progress indicator
            if (mod(trial_ind,25) == 0)
                sprintf('On trail %d',trial_ind)
            end
            
            %% Generate random vector of parameters and then bin them
            
            %per trial channel Parameters (rand draws from uniform 0,1)
            p_act = rand(1,N);
            
                        %Bin each of the actual P's.
            p_act_bin = zeros(1,N);
            for k = 1:1:N
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
            
            %% SPRT INITILIZATION
            
            dbg_print('Start SPRT',DEBUG);
            
            %cost matrix for each channel at each time instant j
            Trial_cost_history = zeros(N,Test_limit);
            
            %path taken by each channel
            Trial_sample_value_hitory = zeros(N,Test_limit);
            
            %distance history
            Trial_dist_history = zeros(N,Test_limit);
            
            %The minimum outer threshold distance
            Trial_e_k = zeros(1,N);
            
            %the current cost
            Trial_cur_cost = zeros(1,N);
            
            %What decision was made for channel K
            Trial_SPRT_decision = zeros(1,N);
            
            %when was it made (what global index)
            Trial_decided_at = zeros(1,N);
            
            %The sequence of Channels that was picked
            Trial_pick_sequence = zeros(1,Test_limit);
            
            %running sum of ones
            Trial_sprt_d_k = zeros(1,N);
            
            %number of measurements each channel k has recieved
            Trial_sprt_m_k = ones(1,N);
            
            %mean number of samples between decision;
            Trail__samples_between_decsion = zeros(1,N);
            
            %Count of samples until a decision is made
            Trial_samples_since_last_decision = 0;
            
            %how many channels below the threshold were touched but not
            %completed.
            Trial_SPRT_mean_channels_passed_on = 0;
            
            %How many samples did those channels that were passed on take?
            Trial_SPRT_mean_channels_passed_on_sample_count = 0;
            
            %track the ammount of time forgivness is callled
            Trial_forgivness_given = zeros(1,N);
            
            %initalize cost
            for k = 1:1:N
                Trial_cur_cost(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,Trial_sprt_d_k(k),Trial_sprt_m_k(k),Test_limit);
            end
            
            %% SPRT Simulation LOOP - The actual SPRT run
            for global_clock = 1:1:Test_limit
                %if All decision have been made stop this run
                if sum(Trial_SPRT_decision > 0 ) == N
                    break;
                end
                
                %bump up the sample count
                Trial_samples_since_last_decision = Trial_samples_since_last_decision + 1;
                
                %pick the channel witht he lowest cost
                cur_pick = pick(Trial_cur_cost);
                
                %record which index I picked
                Trial_pick_sequence(global_clock) = cur_pick;
                
                %sample that channel
                Trial_sprt_d_k(cur_pick) = Trial_sprt_d_k(cur_pick) + random('bino',1,p_act(cur_pick));
                
                %increment it's meausrement count (local clock)
                Trial_sprt_m_k(cur_pick) = Trial_sprt_m_k(cur_pick) + 1;
                
                for k = 1:1:N
                    if  Trial_SPRT_decision(k) == 0
                        if Trial_sprt_d_k(k) > Outer_upper_boundary(Trial_sprt_m_k(k)+1)
                            %decide region 3
                            %            sprintf('Dec:3,chan:%d,sum:%d,bound:%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Outer_upper_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                            Trial_cur_cost(k) = max_cost;
                            Trial_SPRT_decision(k) = 3;
                            %Store the samples between decison count and
                            %reset
                            Trail__samples_between_decsion(k) = Trial_samples_since_last_decision;
                            Trial_samples_since_last_decision = 0;
                            %store the global clock where it was decided
                            Trial_decided_at(k) = global_clock;
                        elseif Trial_sprt_d_k(k) < Outer_lower_boundary(Trial_sprt_m_k(k)+1)
                            %decide region 1
                            %      sprintf('Dec:1,chan:%d,sum:%d,bound:%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Outer_lower_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                            Trial_cur_cost(k) = max_cost;
                            Trial_SPRT_decision(k) = 1;
                            %Store the samples between decison count and
                            %reset
                            Trail__samples_between_decsion(k) = Trial_samples_since_last_decision;
                            Trial_samples_since_last_decision = 0;
                            %store the global clock where it was decided
                            Trial_decided_at(k) = global_clock;
                        elseif Inner_lower_boundary(Trial_sprt_m_k(k)+1) < Trial_sprt_d_k(k) && Trial_sprt_d_k(k) < Inner_upper_boundary(Trial_sprt_m_k(k)+1)
                            %decide region 2
                            %     sprintf('Dec:2,chan:%d,sum:%d,bound:%f_%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Inner_lower_boundary(Trial_sprt_m_k(k)),Inner_upper_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                            Trial_cur_cost(k) = max_cost;
                            Trial_SPRT_decision(k) = 2;
                            %Store the samples between decison count and
                            %reset
                            Trail__samples_between_decsion(k) = Trial_samples_since_last_decision;
                            Trial_samples_since_last_decision = 0;
                            %store the global clock where it was decided
                            Trial_decided_at(k) = global_clock;
                        else
                            %keep sampling
                            %TODO think about this more carefully. This
                            %might be why we never get closer than 0.06
                            if forgive_clock > Trial_sprt_m_k(k)
                                % if we are have only made 2 mistakes above
                                % or below                                
                                if (forgive_allow > Trial_sprt_d_k(k)) || ((Trial_sprt_m_k(k) - Trial_sprt_d_k(k)) < forgive_allow)
                                    forgivness = ForgiveAmt(Outer_upper_boundary,Outer_lower_boundary,Trial_sprt_d_k(k),Trial_sprt_m_k(k),Test_limit);
                                    %Track the ammount of forgiveness used during a trail.
                                    %TODO what the hell is this number
                                    %supposed to be?
                                    Trial_forgivness_given(trial_ind) = Trial_forgivness_given(trial_ind) + 1;
                                end
                            else
                                forgivness = 0;
                            end
                            Trial_cur_cost(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,Trial_sprt_d_k(k),Trial_sprt_m_k(k),Test_limit) + (Trial_sprt_m_k(k) * lambda) + forgivness;
                        end
                    end
                end
                                                
                %Store the trials cost vector (across all channels) on per indexed by time
                Trial_cost_history(:,global_clock) = Trial_cur_cost';
                %Store the trials path vector (across all channels) on per indexed by time
                Trial_sample_value_hitory(:,global_clock) = Trial_sprt_d_k';
                %Store the trials minimum distance vector (across all channels) on per indexed by time
                Trial_dist_history(:,global_clock) = Trial_e_k';
            end
            
            %% SPRT data Analysis - Collecting error and decision information
            
            % Use the bins to compute the error rates
            [Data_SPRT_error_1(lambda_ind,trial_ind,chan_num_ind),Data_SPRT_error_2(lambda_ind,trial_ind,chan_num_ind)] = TrialAlphaBeta(N,p_act_bin,Trial_SPRT_decision);
            dbg_str = sprintf('Trial SPRT error rates are alpha = %f, beta = %f ',Data_SPRT_error_1(lambda_ind,trial_ind,chan_num_ind),Data_SPRT_error_2(lambda_ind,trial_ind,chan_num_ind));
            dbg_print(dbg_str,DEBUG);
            
            %number of channels completed for trial i
            Data_SPRT_Significant(lambda_ind,trial_ind,chan_num_ind) =  sum(Trial_SPRT_decision > 0);
            dbg_str = sprintf('SPRT charaterized %d at deadline, they were:',Data_SPRT_Significant(lambda_ind,trial_ind,chan_num_ind));
            dbg_print(dbg_str,DEBUG);
            
            % Percent Classified by SPRT
            Data_SPRT_fraction_classified(lambda_ind,trial_ind,chan_num_ind) = Data_SPRT_Significant(lambda_ind,trial_ind,chan_num_ind) / N;
            
            Trial_SPRT_p_hat = zeros(1,N);
            %Compute mean p_hat
            for k = 1:1:N
                if Trial_sprt_m_k(k) ~= 0 && Trial_sprt_d_k(k) ~= 0
                    Trial_SPRT_p_hat(k) =  Trial_sprt_d_k(k) / Trial_sprt_m_k(k);
                end
            end
            
            % Percent Average Bin Value
            Data_SPRT_mean_bin_1(lambda_ind,trial_ind,chan_num_ind) = mean(Trial_SPRT_p_hat(Trial_SPRT_decision == 1));
            Data_SPRT_mean_bin_2(lambda_ind,trial_ind,chan_num_ind) = mean(Trial_SPRT_p_hat(Trial_SPRT_decision == 2));
            Data_SPRT_mean_bin_3(lambda_ind,trial_ind,chan_num_ind) = mean(Trial_SPRT_p_hat(Trial_SPRT_decision == 3));
            
            %Mean samples between decsion
            Data_SPRT_mean_samples_between_decsion(lambda_ind,trial_ind,chan_num_ind) = mean(Trail__samples_between_decsion(Trail__samples_between_decsion ~= 0));
            
            %means samples to make a decision
            Data_SPRT_mean_samples_make_decision(lambda_ind,trial_ind,chan_num_ind) = mean(Trial_sprt_m_k(Trial_SPRT_decision ~= 0));
            
            %Mean across all channels (including the zeros) Cuz they don't
            %understand what mean across measured means
            Data_SPRT_mean_m_k(lambda_ind,trial_ind,chan_num_ind) = mean(Trial_sprt_m_k);
            dbg_str = sprintf('Mean number of samples used = %f ',Data_SPRT_mean_m_k(lambda_ind,trial_ind,chan_num_ind));
            dbg_print(dbg_str,DEBUG);
            
            %count the number of channels we've passed on that we shouldn't
            %have
            Data_SPRT_mean_channels_passed_on(lambda_ind,trial_ind,chan_num_ind) = TrialMeanPass(Trial_SPRT_decision,Trial_sprt_m_k,p_act_bin);
            
            %compute the SPRT histogram
            Data_SPRT_histogram(lambda_ind,trial_ind,chan_num_ind,:) = TrialHistogram(Trial_sprt_m_k,Histo_bin);
            
            dbg_print('End SPRT',DEBUG);
            
            %% Sanity Diagonstic output
            if DEBUG
                %         sprintf('Bins')
                % p_act_bin
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
%% Save the results for graphing
matlabpool close;
fname = strcat('DataMultiSprtVaryLlambda-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);

