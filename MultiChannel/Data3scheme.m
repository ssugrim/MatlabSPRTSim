%Generates "Data" for a 3 scheme analysis

%% Paths and local functions

%for use during cmd line run in linux
%addpath '../lib'

%Print diagonstic Messages
DEBUG = true;

%% GLOBAL PARAMETERS  

%Our sample allotment
deadline = 1000;

%number of Trials to test
Trials = 2;

%number of channels 
K = 2^6;

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
Data_good_decision_time = zeros(1,Trials);

%Number of good channels available
Data_good_available = zeros(1,Trials);

%Number of channels that were completed to a decision
Data_comp_full = zeros(1,Trials);

%Errors made during the trial
Data_SPRT_error = zeros(Trials,2);

%Mean number of measurements per channel for SPRT
Data_SPRT_mean_m_k = zeros(1,Trials);

%Errors made with the simple Stratgey
Data_Simple_error = zeros(Trials,2);

%Mean number of measurements per channel for Simple
Data_Simple_mean_m_k = zeros(1,Trials);

%Errors made with the simple Stratgey
Data_Tree_error = zeros(Trials,2);

%Mean number of measurements per channel for TREE
Data_Tree_mean_m_k = zeros(1,Trials);

%% Simulation Begins here 

%%%% Begin Trial Loop (i is my trial index) %%%%
for i = 1:1:Trials 
    
    %% cheap progress indicator
    sprintf('On trail %d',i)
   
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
    Data_good_available(i) = sum(p_act < p_prime_ls);
      
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
    %number of channels completed for trial i
    Data_comp_full(i) =  sum(Trial_SPRT_decision > 0);
    dbg_str = sprintf('SPRT charaterized %d at deadline, they were:',Data_comp_full(i));
    dbg_print(dbg_str,DEBUG);
    
    %Daignostic Display
    %p_act(decision > 0)
    
    %determining all the errors that happened upto the deadline
    
    %Per trial vectors of errors, non zero value means an error occured at that channel index.
    sprt_errors = zeros(1,K);
    
    %Compare Bin Values
    for k = 1:1:K
        if bin_p_act(k) ~= Trial_SPRT_decision(k) && Trial_SPRT_decision(k) ~= 0
            if bin_p_act(k) == 2
                %if the proper bin was 2, and we did not say 2, that is a type 2 error
                sprt_errors(k)=2;
            else
                %other wise it was a type 1 error (False rejection of null %hypothesis)
                sprt_errors(k)=1;
            end
        end
    end
    
    %compute error rates as number of mistakes divided by number classified
    sprt_type1_rate = sum(sprt_errors == 1) / Data_comp_full(i);
    sprt_type2_rate = sum(sprt_errors == 2) / Data_comp_full(i);
    
    %Store and Report error rates
    Data_SPRT_error(i,:) = [sprt_type1_rate,sprt_type2_rate];
    dbg_str = sprintf('SPRT #Type 1 error = %d, #Type 2 error = %d',sum(sprt_errors == 1),sum(sprt_errors == 2));
    dbg_print(dbg_str,DEBUG);
    dbg_str = sprintf('Trial SPRT error rates are alpha = %f, beta = %f ',Data_SPRT_error(i,1),Data_SPRT_error(i,2));
    dbg_print(dbg_str,DEBUG);
    
    %Store and report Mean number of samples used
    Data_SPRT_mean_m_k(i) = mean(Trial_sprt_m_k(Trial_sprt_m_k > 0));
    dbg_str = sprintf('Mean number of samples used = %f ',Data_SPRT_mean_m_k(i));
    dbg_print(dbg_str,DEBUG);
    
    %how long did it take to find 4 good ones.
    Data_good_decision_time(i) = max(Trial_good_decision);
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
    
    %Count the errors made
    simple_errors = zeros(1,K);
    
    %Compare Bin Values
    for k = 1:1:K
        if bin_p_act(k) ~= Trial_simple_bin_hat_p(k) && Trial_simple_m_k(k) ~= 0
%            sprintf('Error Found, hat(p) = %f, p_act = %f, hat_bin = %d, act_bin = %d',Trial_simple_hat_p_k(k), p_act(k), Trial_simple_bin_hat_p(k), bin_p_act(k))
%            sprintf('d_k = %d, m_k = %d, k = %d',Trial_simple_d_k(k),Trial_simple_m_k(k),k)
            if bin_p_act(k) == 2
                %if the proper bin was 2, and we did not say 2, that is a type 2 error
                simple_errors(k)=2;
            else
                %other wise it was a type 1 error (False rejection of null %hypothesis)
                simple_errors(k)=1;
            end
        end
    end
    
    %number of channels with a measurement
    Trial_simple_measured = length(Trial_simple_m_k(Trial_simple_m_k > 0));
    
    dbg_str = sprintf('Number of measured channels = %d', Trial_simple_measured );
    dbg_print(dbg_str,DEBUG);
    %compute error rates as number of mistakes divided by number classified
    simple_type1_rate = sum(simple_errors == 1) / Trial_simple_measured;
    simple_type2_rate = sum(simple_errors == 2) / Trial_simple_measured;
    
    %Store and report Error Rate
    Data_Simple_error(i,:) = [simple_type1_rate,simple_type2_rate];
    dbg_str = sprintf('Simple #Type 1 error = %d, #Type 2 error = %d',sum(simple_errors == 1),sum(simple_errors == 2));
    dbg_print(dbg_str,DEBUG);
    dbg_str = sprintf('Trial Simple error rates are alpha = %f, beta = %f ',Data_Simple_error(i,1),Data_Simple_error(i,2));
    dbg_print(dbg_str,DEBUG);
    
    %Store and report Mean number of samples used
    Data_Simple_mean_m_k(i) = mean(Trial_simple_m_k(Trial_simple_m_k > 0));
    dbg_str = sprintf('Mean number of samples used = %f ',Data_Simple_mean_m_k(i));
    dbg_print(dbg_str,DEBUG);
        
    dbg_print('End Simple',DEBUG);
    
    %% Tree Scheme Initilization
    dbg_print('Start Tree',DEBUG);
    
    %running sum of ones
    Trial_tree_d_k = zeros(1,K);
    
    %number of measurements each channel k has recieved
    Trial_tree_m_k = zeros(1,K);
        
    %Simple scheme estimate of p_k
    Trial_tree_hat_p_k = zeros(1,K);
    
    %bins for the simple scheme
    Trial_tree_bin_hat_p = zeros(1,K);
    
    %slack factor - number of ones we are willing to miss
    slack = 0.25;
    
    %Samples per channel per pass 
    pass_sample = 4;
    
    %% Tree Scheme Simulation Loop
    
    %scan from lowest to highest (taking count of how many measuments each 
    %channel gets). Here we stop at deadline/2 since we're doing two
    %measurements per channel
    
    
    %first pass m is my time index
    if K > deadline
        %if we have more channels than samples, just sample each channel
        %the required ammount for one pass until we run out of samples
        used =  1;
        cur_index = 1;
        while used < deadline + 1
            for n = 1:1:pass_sample
                Trial_tree_d_k(cur_index) = Trial_tree_d_k(cur_index) + random('bino',1,p_act(cur_index));
                Trial_tree_m_k(cur_index) = Trial_tree_m_k(cur_index) + 1;
                used = used + 1;
            end
            cur_index = cur_index + 1;
        end
    else
        %Other wise we do a first pass, and give the required number of
        %samples to each channel.
        
        used =  1;
        for k = 1:1:K
            for n = 1:1:pass_sample
                Trial_tree_d_k(k) = Trial_tree_d_k(k) + random('bino',1,p_act(k));
                Trial_tree_m_k(k) = Trial_tree_m_k(k) + 1;
                used = used + 1;
            end
        end
        
 %       sprintf('Completed first pass, used = %d',used)
        cur_index = 1;
        pass = 0;
        
        %and now distribute measurements that are on the edges
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
                ratio = Trial_tree_d_k(k) / Trial_tree_m_k(k);
                % parameters are within some bound and the channels haven't
                % passed some sane limit
                if (ratio < slack || ratio > 1 - slack) && Trial_tree_m_k(k) < 50
                    picked(k) = 1;
                end
            end

            %Didn't find any thing worth picking, relax the criteria a bit.
            if sum(picked) == 0
                slack = slack + 0.05;
                sprintf('Bumped up the slack to %f', slack)
                continue;
            end
                        
 %           sprintf('Pass = %d, # picked = %d',pass, length(picked(picked > 0)))
            
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
    end
    
    %% Tree Scheme Data Analysis
    
    %computes the estimates of paramters for channels we have measrements for
    for k = 1:1:K
        if Trial_tree_m_k(k) ~= 0
            Trial_tree_hat_p_k(k) = Trial_tree_d_k(k) / Trial_tree_m_k(k);
        end
    end
        
    %bin the Estimates
    for k = 1:1:K
        if Trial_tree_hat_p_k(k) > p_prime_rs
            Trial_tree_bin_hat_p(k) = 3;
        elseif Trial_tree_hat_p_k(k) < p_prime_ls
            Trial_tree_bin_hat_p(k) = 1;
        else
            %if it was never measured, throw it in 2.
            Trial_tree_bin_hat_p(k) = 2;
        end
    end
    
    %Count the errors made
    tree_errors = zeros(1,K);
    
    %Compare Bin Values
    for k = 1:1:K
        if bin_p_act(k) ~= Trial_tree_bin_hat_p(k) && Trial_tree_m_k(k) ~= 0
%            sprintf('Error Found, hat(p) = %f, p_act = %f, hat_bin = %d, act_bin = %d',Trial_simple_hat_p_k(k), p_act(k), Trial_simple_bin_hat_p(k), bin_p_act(k))
%            sprintf('d_k = %d, m_k = %d, k = %d',Trial_simple_d_k(k),Trial_simple_m_k(k),k)
            if bin_p_act(k) == 2
                %if the proper bin was 2, and we did not say 2, that is a type 2 error
                tree_errors(k)=2;
            else
                %other wise it was a type 1 error (False rejection of null %hypothesis)
                tree_errors(k)=1;
            end
        end
    end
  
    
    
    %number of channels with a significant measurement
    Trial_tree_significant_measured = length(Trial_tree_m_k(Trial_tree_m_k > pass_sample));
    dbg_str = sprintf('Number of channels with significant= %d', Trial_tree_significant_measured );
    dbg_print(dbg_str,DEBUG);
    
    Trial_tree_measured = length(Trial_tree_m_k(Trial_tree_m_k > 0));
    dbg_str = sprintf('Number of measured channels = %d', Trial_tree_measured );
    dbg_print(dbg_str,DEBUG);
    
    %compute error rates as number of mistakes divided by number classified
    tree_type1_rate = sum(tree_errors == 1) / Trial_tree_measured;
    tree_type2_rate = sum(tree_errors == 2) / Trial_tree_measured;
    
    %Store and report Error Rate
    Data_Tree_error(i,:) = [tree_type1_rate,tree_type2_rate];
    dbg_str = sprintf('Tree #Type 1 error = %d, #Type 2 error = %d',sum(tree_errors == 1),sum(tree_errors == 2));
    dbg_print(dbg_str,DEBUG);
    dbg_str = sprintf('Trial Tree error rates are alpha = %f, beta = %f ',Data_Tree_error(i,1),Data_Tree_error(i,2));
    dbg_print(dbg_str,DEBUG);
    
    %Store and report Mean number of samples used
    Data_Tree_mean_m_k(i) = mean(Trial_tree_m_k(Trial_tree_m_k > 0));
    dbg_str = sprintf('Mean number of samples used = %f ',Data_Tree_mean_m_k(i));
    dbg_print(dbg_str,DEBUG);
    
    dbg_print('End Tree',DEBUG);
    
    %% Sanity Diagonsitc output
    if DEBUG
%         sprintf('Bins')
%         bin_p_act
%         Trial_SPRT_decision
%         Trial_simple_bin_hat_p
%         Trial_tree_bin_hat_p
%         
%         sprintf('Estimates')
         p_act
%         Trial_sprt_m_k
%         Trial_simple_hat_p_k
%        Trial_simple_d_k
%        Trial_simple_m_k
        Trial_tree_hat_p_k
        Trial_tree_d_k
        Trial_tree_m_k
    end
        
    %% End Trial
end 

%% Store the results for comparison
fname = sprintf('Data3SchemeChan%dT%d.mat',K,Trials);
save(fname);