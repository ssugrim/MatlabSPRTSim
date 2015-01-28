%Generates "Data" for a Multiple parameter (delta, split, lambda) choices of the sprt, evaluating in the same resource regieme L\N 
%Use this to find the optimal 3 paramter choices for a specific resource
%regieme. Can choose various critieara, like max discovered channels, min
%samples to discover, min passed on, etc..

%% Paths and local functions

matlabpool open;

%for use during cmd line run in linux
addpath '../lib'

%Print diagonstic Messages
DEBUG = false;

%% GLOBAL PARAMETERS  

%Our sample allotment
Test_limit = 1000;

%number of Trials to test
Trials = 36;

%Histogram bins
Histo_max = 30;

%Number of Channels to test aganst
Channels = 2^10;

%% SPRT parameter sets

%thresholds parameters
%delta_range = 0.125;
delta_range = 0.05:0.01:0.15
%split_range = 0.3;
split_range = 0.1:0.05:0.9;
p_prime = 0.2;


%alpha and beta design choices
alpha = 0.01;
beta = 0.05;

%penalitly coefficent;
lambda_range =  0:0.05:1;
%lambda_range = [0,0.2,0.5,0.9];

%test line recentering  fudge factor
fudge_factor = 0.98;

%% DATA STORAGE VARS - We do analysis on these


%% SPRT Data

%Errors made during the trial
Data_SPRT_error_1 = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));
Data_SPRT_error_2 = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));

%Mean number of measurements per channel for SPRT
Data_SPRT_mean_m_k = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));

% Fraction Classified by SPRT
Data_SPRT_fraction_classified = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));

% Fraction Classified by SPRT
Data_SPRT_channels_found_H0 = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));
Data_SPRT_channels_found_H1 = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));

%Mean samples between decsion
Data_SPRT_mean_samples_to_decide = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));

%means samples to make a decision
Data_SPRT_mean_samples_make_decision = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));

%means samples to make a decision
Data_SPRT_mean_channels_passed_on = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range));

%SPRT allocation histogram
Data_SPRT_histogram = zeros(length(Channels),length(lambda_range),length(delta_range),length(split_range),(Histo_max + 1));


for lambda_ind  = 1:length(lambda_range)
    for delta_ind = 1:length(delta_range)
        for split_ind = 1:length(split_range)
            %% choice params
            lambda = lambda_range(lambda_ind);
            delta = delta_range(delta_ind);
            split = split_range(split_ind);
            sprintf('On lambda %f,delta %f, split %f',lambda,delta,split)
            
            %% DERVIED VALUES
            %Upper and lower sequential test thresholds (aproximated)
            A = (1 - beta) / alpha;
            B = beta / ( 1 - alpha);
            
            %%%%%%%%%% right side hypothesis %%%%%%%%%%%%%%%
            
            %tuneing parameters - center and width
            p_prime_rs = 1 - p_prime;
            
            %boundaries of accept and reject region
            p_1_rs = p_prime_rs + 2*(delta * (1 - split));
            p_0_rs = p_prime_rs - 2*(delta * split);
            
            %generate test lines
            [L_1_rs,L_0_rs] = GenTestLines(alpha,beta,p_1_rs,p_0_rs,(Test_limit*2));
            
            %%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%
            
            %tuneing parameters - center and width
            p_prime_ls = p_prime;
            
            %boundaries of accept and reject region
            p_1_ls = p_prime_ls + 2*(delta * (1 - split));
            p_0_ls = p_prime_ls - 2*(delta * split);
            
            %generate test lines
            [L_1_ls,L_0_ls] = GenTestLines(alpha,beta,p_1_ls,p_0_ls,(Test_limit*2));
            
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
                
                %Errors made during the trial
                Trial_SPRT_error_1 = zeros(1,Trials);
                Trial_SPRT_error_2 = zeros(1,Trials);
                
                %Mean number of measurements per channel for SPRT
                Trial_SPRT_mean_m_k = zeros(1,Trials);
                
                % Fraction Classified by SPRT
                Trial_SPRT_fraction_classified = zeros(1,Trials);
                
                % Fraction Classified by SPRT
                Trial_SPRT_channels_found_H0 = zeros(1,Trials);
                Trial_SPRT_channels_found_H1 = zeros(1,Trials);
                %Mean samples between decsion
                Trial_SPRT_mean_samples_to_decide = zeros(1,Trials);
                
                %means samples to make a decision
                Trial_SPRT_mean_samples_make_decision = zeros(1,Trials);
                
                %means samples to make a decision
                Trial_SPRT_channels_passed_on = zeros(1,Trials);
                
                %SPRT allocation histogram
                Trial_SPRT_histogram = zeros(Trials,(Histo_max+1));
                
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
                                       
                    %the current cost
                    run_cur_cost = zeros(1,N);
                    
                    %What decision was made for channel K
                    run_SPRT_decision = zeros(1,N);
                    
                    %when was it made (what global index)
                    run_decided_at = zeros(1,N);
                    
                    %The sequence of Channels that was picked
                    run_pick_sequence = zeros(1,Test_limit);
                    
                    %running sum of ones
                    run_sprt_d_k = zeros(1,N);
                    
                    %number of measurements each channel k has recieved
                    run_sprt_m_k = zeros(1,N);
                    
                    %number of spent to decide this channel;
                    run_samples_to_decsion = zeros(1,N);
                    
                    %Count of samples until a decision is made
                    run_samples_since_last_decision = 0;
                    
                    %initalize cost
                    for k = 1:1:N
                        run_cur_cost(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,run_sprt_d_k(k),run_sprt_m_k(k),Test_limit);
                    end
                    
                    %% SPRT Simulation LOOP - The actual SPRT run
                    for global_clock = 1:1:Test_limit
                        %if All decision have been made stop this run
                        if sum(run_SPRT_decision > 0 ) == N
                            break;
                        end
                        
                        %bump up the sample count
                        run_samples_since_last_decision = run_samples_since_last_decision + 1;
                        
                        %pick the channel witht he lowest cost
                        cur_pick = pick(run_cur_cost);
                        
                        %record which index I picked
                        run_pick_sequence(global_clock) = cur_pick;
                        
                        %sample that channel
                        run_sprt_d_k(cur_pick) = run_sprt_d_k(cur_pick) + random('bino',1,p_act(cur_pick));
                        
                        %increment it's meausrement count (local clock)
                        run_sprt_m_k(cur_pick) = run_sprt_m_k(cur_pick) + 1;
                        
                        for k = 1:1:N
                            if  run_SPRT_decision(k) == 0
                                if run_sprt_d_k(k) > Outer_upper_boundary(run_sprt_m_k(k)+1)
                                    %decide region 3
                                    %            sprintf('Dec:3,chan:%d,sum:%d,bound:%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Outer_upper_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                                    run_cur_cost(k) = max_cost;
                                    run_SPRT_decision(k) = 3;
                                    %Store the samples between decison count and
                                    %reset
                                    run_samples_to_decsion(k) = run_samples_since_last_decision;
                                    run_samples_since_last_decision = 0;
                                    %store the global clock where it was decided
                                    run_decided_at(k) = global_clock;
                                elseif run_sprt_d_k(k) < Outer_lower_boundary(run_sprt_m_k(k)+1)
                                    %decide region 1
                                    %      sprintf('Dec:1,chan:%d,sum:%d,bound:%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Outer_lower_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                                    run_cur_cost(k) = max_cost;
                                    run_SPRT_decision(k) = 1;
                                    %Store the samples between decison count and
                                    %reset
                                    run_samples_to_decsion(k) = run_samples_since_last_decision;
                                    run_samples_since_last_decision = 0;
                                    %store the global clock where it was decided
                                    run_decided_at(k) = global_clock;
                                elseif Inner_lower_boundary(run_sprt_m_k(k)+1) < run_sprt_d_k(k) && run_sprt_d_k(k) < Inner_upper_boundary(run_sprt_m_k(k)+1)
                                    %decide region 2
                                    %     sprintf('Dec:2,chan:%d,sum:%d,bound:%f_%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Inner_lower_boundary(Trial_sprt_m_k(k)),Inner_upper_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                                    run_cur_cost(k) = max_cost;
                                    run_SPRT_decision(k) = 2;
                                    %Store the samples between decison count and
                                    %reset
                                    run_samples_to_decsion(k) = run_samples_since_last_decision;
                                    run_samples_since_last_decision = 0;
                                    %store the global clock where it was decided
                                    run_decided_at(k) = global_clock;
                                else
                                    %keep sampling
                                   run_cur_cost(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,run_sprt_d_k(k),run_sprt_m_k(k),Test_limit) + (run_sprt_m_k(k) * lambda);
                                end
                            end
                        end
                    end                   
                    
                    %Compute type 1 and type2 error
                    [Trial_SPRT_error_1(trial_ind),Trial_SPRT_error_2(trial_ind)] = TrialAlphaBeta(N,p_act_bin,run_SPRT_decision);
                    dbg_str = sprintf('Trial SPRT error rates are alpha = %f, beta = %f ',Trial_SPRT_error_1(trial_ind),Trial_SPRT_error_2(trial_ind));
                    dbg_print(dbg_str,DEBUG);
                    
                    %channels found
                    Trial_SPRT_channels_found_H0(trial_ind) = sum(run_SPRT_decision == 1) + sum(run_SPRT_decision == 3);
                    Trial_SPRT_channels_found_H1(trial_ind) = sum(run_SPRT_decision == 2);
        
                    % Percent Classified by SPRT
                    Trial_SPRT_fraction_classified(trial_ind) = sum(run_SPRT_decision ~= 0) / N;
                                      
                    %Mean samples between decsion
                    Trial_SPRT_mean_samples_to_decide(trial_ind) = mean(run_samples_to_decsion(run_samples_to_decsion ~= 0));
                    
                    %means samples to make a decision
                    Trial_SPRT_mean_samples_make_decision(trial_ind) = mean(run_sprt_m_k(run_SPRT_decision ~= 0));
                    
                    %Mean across all channels which were tested.
                    Trial_SPRT_mean_m_k(trial_ind) = mean(run_sprt_m_k(run_sprt_m_k ~= 0));
                    dbg_str = sprintf('Mean number of samples used = %f ',Trial_SPRT_mean_m_k(trial_ind));
                    dbg_print(dbg_str,DEBUG);
                    
                    %count the number of channels we've passed on that we shouldn't have
                    Trial_SPRT_channels_passed_on(trial_ind) = TrialChanPass(run_SPRT_decision,run_sprt_m_k,p_act_bin);
                    
                    %compute the SPRT histogram
                    Histo_bin = zeros(1,Histo_max+1);
                    for hist_ind = 0:1:Histo_max
                        Histo_bin(hist_ind+1) = sum(run_sprt_m_k == hist_ind);
                    end
                    
                    Trial_SPRT_histogram(trial_ind,:) = Histo_bin;
                    
                    dbg_print('End SPRT',DEBUG);
                end
                
                %% SPRT data Analysis - Collecting error and decision information
                % Use the bins to compute the error rates
                %Errors made during the trial
                Data_SPRT_error_1(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_error_1);
                Data_SPRT_error_2(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_error_2);
                
                %Mean number of measurements per channel for SPRT
                Data_SPRT_mean_m_k(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_mean_m_k);
                
                % Fraction Classified by SPRT
                Data_SPRT_fraction_classified(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_fraction_classified);
                
                % Fraction Classified by SPRT
                Data_SPRT_channels_found_H0(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_channels_found_H0);
                Data_SPRT_channels_found_H1(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_channels_found_H1);
                
                %Mean samples between decsion
                Data_SPRT_mean_samples_to_decide(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_mean_samples_to_decide);
                
                %means samples to make a decision
                Data_SPRT_mean_samples_make_decision(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_mean_samples_make_decision);
                
                %means samples to make a decision
                Data_SPRT_mean_channels_passed_on(chan_num_ind,lambda_ind,delta_ind,split_ind) = mean(Trial_SPRT_channels_passed_on);
                
                %SPRT allocation histogram
                Mean_Histo = zeros(1,Histo_max+1);
                %compute the mean across trials for each bin.
                for hist_ind = 1:1:Histo_max+1
                    Mean_Histo(hist_ind) = mean(Trial_SPRT_histogram(:,hist_ind));
                end                
                Data_SPRT_histogram(chan_num_ind,lambda_ind,delta_ind,split_ind,:) = Mean_Histo;
            end
        end
    end
end

%% Save the results for graphing
matlabpool close;
fname = strcat('DataMultiSprtVaryDeltaSplitLlambda-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);

