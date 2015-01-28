%Generates "Data" for a Multiple parameter (delta, split, lambda) choices of the sprt, evaluating in the same resource regieme L\N but with smaller N for speed. 
%Will have to test in the bigger case later

%% Paths and local functions

matlabpool open;

%for use during cmd line run in linux
%addpath '../lib'

%Print diagonstic Messages
DEBUG = false;

%% GLOBAL PARAMETERS  

%Our sample allotment
Test_limit = 1000;

%number of Trials to test
Trials = 100;

%Histogram bins
Histo_max = 30;

%Number of Channels to test aganst
Channels = 100:100:1000;

%% SPRT parameter sets

%thresholds parameters
delta = 0.15;
%delta_range = 0.125;
split = 0.1;
%split_range = 0.5;
p_prime = 0.2;


%alpha and beta design choices
alpha = 0.01;
beta = 0.05;

%penalitly coefficent;
lambda =  0.3;

%test line recentering  fudge factor
fudge_factor = 0.98;

%% Tree Parameters

 %fraction of samples to keep in reserve for passes after first
reserve = 0.5;

pass_sample = 4;

slack_init = 1 / pass_sample;

%% DATA STORAGE VARS - We do analysis on these


%% SPRT Data

%Errors made during the trial
Data_SPRT_error_1 = zeros(1,length(Channels));
Data_SPRT_error_2 = zeros(1,length(Channels));

%Mean number of measurements per channel for SPRT
Data_SPRT_mean_m_k = zeros(1,length(Channels));

% Fraction Classified by SPRT
Data_SPRT_fraction_classified = zeros(1,length(Channels));

% Fraction Classified by SPRT
Data_SPRT_channels_found_H0 = zeros(1,length(Channels));
Data_SPRT_channels_found_H1 = zeros(1,length(Channels));

%Mean samples between decsion
Data_SPRT_mean_samples_to_decide = zeros(1,length(Channels));

%means samples to make a decision
Data_SPRT_mean_samples_make_decision = zeros(1,length(Channels));

%means samples to make a decision
Data_SPRT_mean_channels_passed_on = zeros(1,length(Channels));

%mean time to complete as number of channels increase
Data_SPRT_CPU_time = zeros(1,length(Channels));

%SPRT allocation histogram
Data_SPRT_histogram = zeros(length(Channels),(Histo_max + 1));


%% Tree Data 
%Errors made during the trial
Data_Tree_error_1 = zeros(1,length(Channels));
Data_Tree_error_2 = zeros(1,length(Channels));

%Mean number of measurements per channel for SPRT
Data_Tree_mean_m_k = zeros(1,length(Channels));

% Fraction Classified by SPRT
Data_Tree_fraction_classified = zeros(1,length(Channels));

% Fraction Classified by SPRT
Data_Tree_channels_found_H0 = zeros(1,length(Channels));
Data_Tree_channels_found_H1 = zeros(1,length(Channels));

%Mean samples between decsion
Data_Tree_mean_samples_to_decide = zeros(1,length(Channels));

%means samples to make a decision
Data_Tree_mean_samples_make_decision = zeros(1,length(Channels));

%mean time to complete as number of channels increase
Data_Tree_CPU_time = zeros(1,length(Channels));

%Tree allocation histogram
Data_Tree_histogram = zeros(length(Channels),(Histo_max + 1));

%% Simple Data 
%Errors made during the trial
Data_Simple_error_1 = zeros(1,length(Channels));
Data_Simple_error_2 = zeros(1,length(Channels));

%Mean number of measurements per channel for SPRT
Data_Simple_mean_m_k = zeros(1,length(Channels));

% Fraction Classified by SPRT
Data_Simple_fraction_classified = zeros(1,length(Channels));

% Fraction Classified by SPRT
Data_Simple_channels_found_H0 = zeros(1,length(Channels));
Data_Simple_channels_found_H1 = zeros(1,length(Channels));

%Mean samples between decsion
Data_Simple_mean_samples_to_decide = zeros(1,length(Channels));

%means samples to make a decision
Data_Simple_mean_samples_make_decision = zeros(1,length(Channels));

%mean time to complete as number of channels increase
Data_Simple_CPU_time = zeros(1,length(Channels));

%Simple allocation histogram
Data_Simple_histogram = zeros(length(Channels),(Histo_max + 1));

           
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
[L_1_rs,L_0_rs] = GenTestLines(alpha,beta,p_1_rs,p_0_rs,(Test_limit*2));

%%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_ls = p_prime;

%boundaries of accept and reject region
p_1_ls = p_prime_ls + (delta * (1 - split));
p_0_ls = p_prime_ls - (delta * split);

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

%% Channel Loop
for chan_num_ind = 1:length(Channels)
    N = Channels(chan_num_ind);
    sprintf('On channel set %d',N)
    
    %% SPRT Trial Data for channel SET N
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
    
    %Trial clock of run time
    Trial_SPRT_CPU_time = zeros(1,Trials);
    
    %SPRT allocation histogram
    Trial_SPRT_histogram = zeros(Trials,(Histo_max+1));
    
    %% Tree Trial Data for channel SET N
    %Errors made during the trial
    Trial_Tree_error_1 = zeros(1,Trials);
    Trial_Tree_error_2 = zeros(1,Trials);
    
    %Mean number of measurements per channel for SPRT
    Trial_Tree_mean_m_k = zeros(1,Trials);
    
    % Fraction Classified by SPRT
    Trial_Tree_fraction_classified = zeros(1,Trials);
    
    % Fraction Classified by SPRT
    Trial_Tree_channels_found_H0 = zeros(1,Trials);
    Trial_Tree_channels_found_H1 = zeros(1,Trials);
    %Mean samples between decsion
    Trial_Tree_mean_samples_to_decide = zeros(1,Trials);
    
    %means samples to make a decision
    Trial_Tree_mean_samples_make_decision = zeros(1,Trials);
    
    %Trial clock of run time
    Trial_Tree_CPU_time = zeros(1,Trials);
    
    %Tree allocation histogram
    Trial_Tree_histogram = zeros(Trials,(Histo_max+1));
    
    %% Simple Trial Data for channel SET N
    %Errors made during the trial
    Trial_Simple_error_1 = zeros(1,Trials);
    Trial_Simple_error_2 = zeros(1,Trials);
    
    %Mean number of measurements per channel for SPRT
    Trial_Simple_mean_m_k = zeros(1,Trials);
    
    % Fraction Classified by SPRT
    Trial_Simple_fraction_classified = zeros(1,Trials);
    
    % Fraction Classified by SPRT
    Trial_Simple_channels_found_H0 = zeros(1,Trials);
    Trial_Simple_channels_found_H1 = zeros(1,Trials);
    %Mean samples between decsion
    Trial_Simple_mean_samples_to_decide = zeros(1,Trials);
    
    %means samples to make a decision
    Trial_Simple_mean_samples_make_decision = zeros(1,Trials);
    
    %Trial clock of run time
    Trial_Simple_CPU_time = zeros(1,Trials);
    
    %Simple allocation histogram
    Trial_Simple_histogram = zeros(Trials,(Histo_max+1));
    
    %% Begin Trial Loop
    parfor trial_ind = 1:Trials
  %  for trial_ind = 1:Trials
        %% cheap progress indicator
        if (mod(trial_ind,25) == 0)
            sprintf('On trail %d',trial_ind)
        end
        
        %% Generate random vector of parameters and then bin them
        
        %per trial channel Parameters (rand draws from uniform 0,1)
        run_p_act = rand(1,N);
        
        run_samples = zeros(N,Test_limit);
        
        %Each classifier gets the same samples and makes different decisons
        %on them (so we're compareing the same thing).
        for k = 1:1:N
            run_samples(k,:) = random('bino',1,run_p_act(k),[1,Test_limit]);
        end
        
        %Bin each of the actual P's.
        run_p_act_bin = zeros(1,N);
        for k = 1:1:N
            if run_p_act(k) > p_prime_rs
                %above the upper limit, must belong to 3
                run_p_act_bin(k) = 3;
            elseif run_p_act(k) < p_prime_ls
                %below the lower limit must belong to 1
                run_p_act_bin(k) = 1;
            else
                %every thing else is in 2
                run_p_act_bin(k) = 2;
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
        run_sprt_samples_to_decsion = zeros(1,N);
        
        %Count of samples until a decision is made
        run_sprt_samples_since_last_decision = 0;
        
        %initalize cost
        for k = 1:1:N
            run_cur_cost(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,run_sprt_d_k(k),run_sprt_m_k(k),Test_limit);
        end
        
        %% SPRT Simulation LOOP - The actual SPRT run
        run_sprt_start_time = cputime;
        
        for global_clock = 1:1:Test_limit
            %if All decision have been made stop this run
            if sum(run_SPRT_decision > 0 ) == N
                break;
            end
            
            %bump up the sample count
            run_sprt_samples_since_last_decision = run_sprt_samples_since_last_decision + 1;
            
            %pick the channel witht he lowest cost
            cur_pick = pick(run_cur_cost);
            
            %record which index I picked
            run_pick_sequence(global_clock) = cur_pick;
            
            %sample that channel, pull samples from the
            run_sprt_d_k(cur_pick) = run_sprt_d_k(cur_pick) + run_samples(cur_pick,run_sprt_m_k(cur_pick)+1);
            
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
                        run_sprt_samples_to_decsion(k) = run_sprt_samples_since_last_decision;
                        run_sprt_samples_since_last_decision = 0;
                        %store the global clock where it was decided
                        run_decided_at(k) = global_clock;
                    elseif run_sprt_d_k(k) < Outer_lower_boundary(run_sprt_m_k(k)+1)
                        %decide region 1
                        %      sprintf('Dec:1,chan:%d,sum:%d,bound:%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Outer_lower_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                        run_cur_cost(k) = max_cost;
                        run_SPRT_decision(k) = 1;
                        %Store the samples between decison count and
                        %reset
                        run_sprt_samples_to_decsion(k) = run_sprt_samples_since_last_decision;
                        run_sprt_samples_since_last_decision = 0;
                        %store the global clock where it was decided
                        run_decided_at(k) = global_clock;
                    elseif Inner_lower_boundary(run_sprt_m_k(k)+1) < run_sprt_d_k(k) && run_sprt_d_k(k) < Inner_upper_boundary(run_sprt_m_k(k)+1)
                        %decide region 2
                        %     sprintf('Dec:2,chan:%d,sum:%d,bound:%f_%f,act:%f,mea:%d,cost:%f',k,Trial_sprt_d_k(k),Inner_lower_boundary(Trial_sprt_m_k(k)),Inner_upper_boundary(Trial_sprt_m_k(k)),p_act(k),Trial_sprt_m_k(k),costs(k))
                        run_cur_cost(k) = max_cost;
                        run_SPRT_decision(k) = 2;
                        %Store the samples between decison count and
                        %reset
                        run_sprt_samples_to_decsion(k) = run_sprt_samples_since_last_decision;
                        run_sprt_samples_since_last_decision = 0;
                        %store the global clock where it was decided
                        run_decided_at(k) = global_clock;
                    else
                        %keep sampling
                        run_cur_cost(k) = StepsToComplete(Outer_upper_boundary,Outer_lower_boundary,run_sprt_d_k(k),run_sprt_m_k(k),Test_limit) + (run_sprt_m_k(k) * lambda);
                    end
                end
            end
        end
        run_sprt_stop_time = cputime;
        
        %% SPRT data Analysis
        %Compute type 1 and type2 error
        [Trial_SPRT_error_1(trial_ind),Trial_SPRT_error_2(trial_ind)] = TrialAlphaBeta(N,run_p_act_bin,run_SPRT_decision);
        dbg_str = sprintf('Trial SPRT error rates are alpha = %f, beta = %f ',Trial_SPRT_error_1(trial_ind),Trial_SPRT_error_2(trial_ind));
        dbg_print(dbg_str,DEBUG);
        
        %channels found
        Trial_SPRT_channels_found_H0(trial_ind) = sum(run_SPRT_decision == 1) + sum(run_SPRT_decision == 3);
        Trial_SPRT_channels_found_H1(trial_ind) = sum(run_SPRT_decision == 2);
        
        % Percent Classified by SPRT
        Trial_SPRT_fraction_classified(trial_ind) = sum(run_SPRT_decision ~= 0) / N;
        
        %Mean samples between decsion
        Trial_SPRT_mean_samples_to_decide(trial_ind) = mean(run_sprt_samples_to_decsion(run_sprt_samples_to_decsion ~= 0));
        
        %means samples to make a decision
        Trial_SPRT_mean_samples_make_decision(trial_ind) = mean(run_sprt_m_k(run_SPRT_decision ~= 0));
        
        %Mean across all channels which were tested.
        Trial_SPRT_mean_m_k(trial_ind) = mean(run_sprt_m_k(run_sprt_m_k ~= 0));
        dbg_str = sprintf('Mean number of samples used = %f ',Trial_SPRT_mean_m_k(trial_ind));
        dbg_print(dbg_str,DEBUG);
        
        %count the number of channels we've passed on that we shouldn't have
        Trial_SPRT_channels_passed_on(trial_ind) = TrialChanPass(run_SPRT_decision,run_sprt_m_k,run_p_act_bin);
        
        %How long did the categorization take?
        Trial_SPRT_CPU_time(trial_ind) = run_sprt_stop_time - run_sprt_start_time;
        
        %compute the SPRT histogram
        Histo_bin = zeros(1,Histo_max+1);
        for hist_ind = 0:1:Histo_max
            Histo_bin(hist_ind+1) = sum(run_sprt_m_k == hist_ind);
        end
        
        Trial_SPRT_histogram(trial_ind,:) = Histo_bin;
        
        
        dbg_print('End SPRT',DEBUG);
        
        %% Tree initilization
        dbg_print('Start Tree',DEBUG);
        
        %must reset the slack each trial.
        slack = slack_init;
        
        %running sum of ones
        run_Tree_d_k = zeros(1,N);
        
        %number of measurements each channel k has recieved
        run_Tree_m_k = zeros(1,N);
        
        %% Tree Simulation Loop
        
        %scan from lowest to highest (taking count of how many measuments each
        %channel gets). Here we stop at Test_limit/2 since we're doing two
        %measurements per channel
        
        %do a first pass with ((1 - reserve) * limit) channels alloted
        
        %count of samples used
        used =  0;
        %first pass channel allotment
        alloted_channels = floor(((1 - reserve) * Test_limit) / pass_sample);
        
        for k = 1:1:min(alloted_channels,N)
            for n = 1:1:pass_sample
                %pull samples from the stored sample matrix
                run_Tree_d_k(k) = run_Tree_d_k(k) + run_samples(k,run_Tree_m_k(k)+1);
                run_Tree_m_k(k) = run_Tree_m_k(k) + 1;
                used = used + 1;
            end
        end
               
        %and now distribute measurements that are on the edges
        pass = 0;
        
        run_Tree_start_time = cputime;        
        while used < Test_limit + 1
            %relaxed too much I must be done (have to let it get past .5
            %otherwise we miss some channels even though we have
            %measurements
            if slack > 0.6
                break;
            end
            
            picked = zeros(1,N);
            %pick a set of channels
            for k = 1:1:N
                if run_Tree_m_k(k) > 0
                    ratio = run_Tree_d_k(k) / run_Tree_m_k(k);
                    % parameters are within some bound and the channels haven't
                    % passed some sane limit
                    if (ratio < slack || ratio > 1 - slack) && run_Tree_m_k(k) < 50
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
            for k = 1:1:N
                if picked(k) ~= 0
                    for n = 1:1:pass_sample
                        run_Tree_d_k(k) = run_Tree_d_k(k) + run_samples(k,run_Tree_m_k(k)+1);
                        run_Tree_m_k(k) = run_Tree_m_k(k) + 1;
                        used = used + 1;
                    end
                end
            end
            
            %pass count
            pass = pass + 1;
        end
        
        
        %Simple scheme estimate of p_k
        run_Tree_hat_p_k = zeros(1,N);
        
        %computes the estimates of paramters for channels we have measrements for
        for k = 1:1:N
            if run_Tree_m_k(k) > 0
                run_Tree_hat_p_k(k) = run_Tree_d_k(k) / run_Tree_m_k(k);
            else
                run_Tree_hat_p_k(k) = -1;
            end
        end
        
        %bins for the simple scheme
        run_Tree_decision = zeros(1,N);
        
        %bin the Estimates
        for k = 1:1:N
            if run_Tree_hat_p_k(k) ~= -1
                if run_Tree_hat_p_k(k) > p_prime_rs
                    run_Tree_decision(k) = 3;
                elseif run_Tree_hat_p_k(k) < p_prime_ls
                    run_Tree_decision(k) = 1;
                else
                    %if it was never measured, throw it in 2.
                    run_Tree_decision(k) = 2;
                end
            end
        end
        run_Tree_stop_time = cputime;
        
        %% Tree Scheme Data Analysis
        [Trial_Tree_error_1(trial_ind),Trial_Tree_error_2(trial_ind)] = TrialAlphaBeta(N,run_p_act_bin,run_Tree_decision);
        dbg_str = sprintf('Trial Tree error rates are alpha = %f, beta = %f',Trial_Tree_error_1(trial_ind),Trial_Tree_error_2(trial_ind));
        dbg_print(dbg_str,DEBUG);
        
        %channels found
        Trial_Tree_channels_found_H0(trial_ind) = sum(run_Tree_decision == 1) + sum(run_Tree_decision == 3);
        Trial_Tree_channels_found_H1(trial_ind) = sum(run_Tree_decision == 2);
        
        % Percent Classified by Tree
        Trial_Tree_fraction_classified(trial_ind) = sum(run_Tree_decision ~= 0) / N;
        
        %Mean samples between decsion (in this case just what ever made it
        %past all first pass. (Perhaps this should be pass all phases?)
        Trial_Tree_mean_samples_to_decide(trial_ind) = mean(run_Tree_m_k(run_Tree_m_k > pass_sample));
        
        %means samples to make a decision
        Trial_Tree_mean_samples_make_decision(trial_ind) = mean(run_Tree_m_k(run_Tree_decision ~= 0));
        
        %Mean across all channels which were tested.
        Trial_Tree_mean_m_k(trial_ind) = mean(run_Tree_m_k(run_Tree_m_k ~= 0));
        dbg_str = sprintf('Mean number of samples used = %f ',Trial_Tree_mean_m_k(trial_ind));
        dbg_print(dbg_str,DEBUG);
        
        %How long did the categorization take?
        Trial_Tree_CPU_time(trial_ind) = run_Tree_stop_time - run_Tree_start_time;
        
        %compute the Tree histogram
        Histo_bin = zeros(1,Histo_max+1);
        for hist_ind = 0:1:Histo_max
            Histo_bin(hist_ind+1) = sum(run_Tree_m_k == hist_ind);
        end
                
        Trial_Tree_histogram(trial_ind,:) = Histo_bin;
        
        dbg_print('End Tree',DEBUG);
        
        %% Simple Initilization
        dbg_print('Start Simple',DEBUG);
        
        %running sum of ones
        run_Simple_d_k = zeros(1,N);
        
        %number of measurements each channel k has recieved
        run_Simple_m_k = zeros(1,N);
        
        %% Simple Scheme Simulation Loop
        
        run_Simple_start_time = cputime;
        
        %Simple scheme estimate of p_k
        run_Simple_hat_p_k = zeros(1,N);
        
        %bins for the simple scheme
        run_Simple_decision = zeros(1,N);
        
        %scan from lowest to highest (taking count of how many measuments each
        %channel gets.
        cur_index = 1;
        for l = 1:1:Test_limit+1
            %Index wrap around
            cur_index = mod(l,N) + 1;
            
            %take a single sample
            run_Simple_d_k(cur_index) = run_Simple_d_k(cur_index) + run_samples(cur_index,run_Simple_m_k(cur_index)+1);
            
            %Count how many times this channel has been measured
            run_Simple_m_k(cur_index) = run_Simple_m_k(cur_index) + 1;
        end
        
        %computes the estimates we have measrements for
        for k = 1:1:N
            if run_Simple_m_k(k) ~= 0
                run_Simple_hat_p_k(k) = run_Simple_d_k(k) / run_Simple_m_k(k);
            end
        end
        
        %bin the Estimates
        for k = 1:1:N
            if run_Simple_m_k(k) ~= 0
                if run_Simple_hat_p_k(k) > p_prime_rs
                    run_Simple_decision(k) = 3;
                elseif run_Simple_hat_p_k(k) < p_prime_ls
                    run_Simple_decision(k) = 1;
                else
                    %if it was never measured, throw it in 2.
                    run_Simple_decision(k) = 2;
                end
            end
        end
        
        run_Simple_stop_time = cputime;
        
        %% Simple Scheme Data Analysis
        
        [Trial_Simple_error_1(trial_ind),Trial_Simple_error_2(trial_ind)] = TrialAlphaBeta(N,run_p_act_bin,run_Simple_decision);
        dbg_str = sprintf('Trial Simple error rates are alpha = %f, beta = %f',Trial_Simple_error_1(trial_ind),Trial_Simple_error_2(trial_ind));
        dbg_print(dbg_str,DEBUG);
        
        %channels found
        Trial_Simple_channels_found_H0(trial_ind) = sum(run_Simple_decision == 1) + sum(run_Simple_decision == 3);
        Trial_Simple_channels_found_H1(trial_ind) = sum(run_Simple_decision == 2);
        
        % Percent Classified by Simple
        Trial_Simple_fraction_classified(trial_ind) = sum(run_Simple_decision ~= 0) / N;
        
        %Mean samples between decsion (in this case just what ever made it
        %past all first pass. (Perhaps this should be pass all phases?)
        Trial_Simple_mean_samples_to_decide(trial_ind) = mean(run_Simple_m_k(run_Simple_m_k > (Test_limit / N)));
        
        %means samples to make a decision
        Trial_Simple_mean_samples_make_decision(trial_ind) = mean(run_Simple_m_k(run_Simple_decision ~= 0));
        
        %Mean across all channels which were tested.
        Trial_Simple_mean_m_k(trial_ind) = mean(run_Simple_m_k(run_Simple_m_k ~= 0));
        dbg_str = sprintf('Mean number of samples used = %f ',Trial_Simple_mean_m_k(trial_ind));
        dbg_print(dbg_str,DEBUG);
        
        %How long did the categorization take?
        Trial_Simple_CPU_time(trial_ind) = run_Simple_stop_time - run_Simple_start_time;
        
        %compute the Simple histogram
        Histo_bin = zeros(1,Histo_max+1);
        for hist_ind = 0:1:Histo_max
            Histo_bin(hist_ind+1) = sum(run_Simple_m_k == hist_ind);
        end
                
        Trial_Simple_histogram(trial_ind,:) = Histo_bin;
        
        dbg_print('End Tree',DEBUG);
        
        
    end
    
    %% SPRT Data Storage - Record means per Channel
    % Use the bins to compute the error rates
    %Errors made during the trial
    Data_SPRT_error_1(chan_num_ind) = mean(Trial_SPRT_error_1);
    Data_SPRT_error_2(chan_num_ind) = mean(Trial_SPRT_error_2);
    
    %Mean number of measurements per channel for SPRT
    Data_SPRT_mean_m_k(chan_num_ind) = mean(Trial_SPRT_mean_m_k);
    
    % Fraction Classified by SPRT
    Data_SPRT_fraction_classified(chan_num_ind) = mean(Trial_SPRT_fraction_classified);
    
    % Fraction Classified by SPRT
    Data_SPRT_channels_found_H0(chan_num_ind) = mean(Trial_SPRT_channels_found_H0);
    Data_SPRT_channels_found_H1(chan_num_ind) = mean(Trial_SPRT_channels_found_H1);
    
    %Mean samples between decsion
    Data_SPRT_mean_samples_to_decide(chan_num_ind) = mean(Trial_SPRT_mean_samples_to_decide);
    
    %means samples to make a decision
    Data_SPRT_mean_samples_make_decision(chan_num_ind) = mean(Trial_SPRT_mean_samples_make_decision);
    
    %means samples to make a decision
    Data_SPRT_mean_channels_passed_on(chan_num_ind) = mean(Trial_SPRT_channels_passed_on);
    
    %means CPU time    
    Data_SPRT_CPU_time(chan_num_ind) = mean(Trial_SPRT_CPU_time);
    
    %SPRT allocation histogram
    Mean_Histo = zeros(1,Histo_max+1);
    
    %compute the mean across trials for each bin.
    for hist_ind = 1:1:Histo_max+1
        Mean_Histo(hist_ind) = mean(Trial_SPRT_histogram(:,hist_ind));
    end
    
    Data_SPRT_histogram(chan_num_ind,:) = Mean_Histo;
    
    %% Tree Data Storage - Record means per Channel
    % Use the bins to compute the error rates
    %Errors made during the trial
    Data_Tree_error_1(chan_num_ind) = mean(Trial_Tree_error_1);
    Data_Tree_error_2(chan_num_ind) = mean(Trial_Tree_error_2);
    
    %Mean number of measurements per channel for Tree
    Data_Tree_mean_m_k(chan_num_ind) = mean(Trial_Tree_mean_m_k);
    
    % Fraction Classified by Tree
    Data_Tree_fraction_classified(chan_num_ind) = mean(Trial_Tree_fraction_classified);
    
    % Fraction Classified by Tree
    Data_Tree_channels_found_H0(chan_num_ind) = mean(Trial_Tree_channels_found_H0);
    Data_Tree_channels_found_H1(chan_num_ind) = mean(Trial_Tree_channels_found_H1);
    
    %Mean samples between decsion
    Data_Tree_mean_samples_to_decide(chan_num_ind) = mean(Trial_Tree_mean_samples_to_decide);
    
    %means samples to make a decision
    Data_Tree_mean_samples_make_decision(chan_num_ind) = mean(Trial_Tree_mean_samples_make_decision);
    
    %means CPU time
    Data_Tree_CPU_time(chan_num_ind) = mean(Trial_Tree_CPU_time);
            
    %Tree allocation histogram
    Mean_Histo = zeros(1,Histo_max+1);
    %compute the mean across trials for each bin.
    for hist_ind = 1:1:Histo_max+1
        Mean_Histo(hist_ind) = mean(Trial_Tree_histogram(:,hist_ind));
    end
    Data_Tree_histogram(chan_num_ind,:) = Mean_Histo;
    
     %% Simple Data Storage - Record means per Channel
    % Use the bins to compute the error rates
    %Errors made during the trial
    Data_Simple_error_1(chan_num_ind) = mean(Trial_Simple_error_1);
    Data_Simple_error_2(chan_num_ind) = mean(Trial_Simple_error_2);
    
    %Mean number of measurements per channel for Simple
    Data_Simple_mean_m_k(chan_num_ind) = mean(Trial_Simple_mean_m_k);
    
    % Fraction Classified by Simple
    Data_Simple_fraction_classified(chan_num_ind) = mean(Trial_Simple_fraction_classified);
    
    % Fraction Classified by Simple
    Data_Simple_channels_found_H0(chan_num_ind) = mean(Trial_Simple_channels_found_H0);
    Data_Simple_channels_found_H1(chan_num_ind) = mean(Trial_Simple_channels_found_H1);
    
    %Mean samples between decsion
    Data_Simple_mean_samples_to_decide(chan_num_ind) = mean(Trial_Simple_mean_samples_to_decide);
    
    %means samples to make a decision
    Data_Simple_mean_samples_make_decision(chan_num_ind) = mean(Trial_Simple_mean_samples_make_decision);
    
    %means CPU time
    Data_Simple_CPU_time(chan_num_ind) = mean(Trial_Simple_CPU_time);
            
    %Simple allocation histogram
    Mean_Histo = zeros(1,Histo_max+1);
    %compute the mean across trials for each bin.
    for hist_ind = 1:1:Histo_max+1
        Mean_Histo(hist_ind) = mean(Trial_Simple_histogram(:,hist_ind));
    end
    Data_Simple_histogram(chan_num_ind,:) = Mean_Histo;
    
end

%% Save the results for graphing
matlabpool close;
fname = strcat('Data3SchemeNewDist-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);