%Generates "Data" for a K channel reward model. 

%for use during cmd line run in linux
%addpath '../lib'

%global Parameters
%number of Trials
Trials = 1;

%number of channels 
K = 8;

%number of good channels we want to dig out
good = 4;

%dealine to count Completed (we'll let the test run to completeion, we
%merely count here).
deadline = 1000;

%penalitly coefficent;
lambda = 0.05;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

delta = 0.1;
p_prime = 0.2;
split = 0.3;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%The hard limt All tests must stop at (Lenght of our vectors, should be
%longer than eta, for some allowance).
Test_Limit = 100000;

%%%%%%%%%% right side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_rs = 1 - p_prime;

%boundaries of accept and reject region
p_1_rs = p_prime_rs + (delta * (1 - split));
p_0_rs = p_prime_rs - (delta * split);

%generate test lines
[L_1_rs,L_0_rs] = GenTestLines(alpha,beta,p_1_rs,p_0_rs,Test_Limit);

%%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_ls = p_prime;

%boundaries of accept and reject region
p_1_ls = p_prime_ls + (delta * (1 - split));
p_0_ls = p_prime_ls - (delta * split);

%generate test lines
[L_1_ls,L_0_ls] = GenTestLines(alpha,beta,p_1_ls,p_0_ls,Test_Limit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convinence Names

%Upper threshold
T_U = L_1_rs;

%Lower threshold
T_L = L_0_ls;

%Center Upper threshold
T_C_U = L_0_rs;

%Center Lower threshold
T_C_L = L_1_ls;

%Trails storage vars
full_decision_time = zeros(K,Trials);
good_decision_time = zeros(1,Trials);
good_available = zeros(1,Trials);
deadline_comp_full = zeros(1,Trials);
Trial_error = zeros(Trials,2);
Dead_Trial_error = zeros(Trials,2);
max_cost =  pdist([Test_Limit,T_U(Test_Limit);Test_Limit,T_L(Test_Limit)],'euclidean') + (lambda * Test_Limit); 

%Begin Simulation

for i = 1:1:Trials
    %cheap progress indicator
    sprintf('On trail %d',i)
    
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
    good_available(i) = sum(p_act < p_prime_ls);
    
    %cost matrix for each channel for at each time instant j
    cost = zeros(K,Test_Limit);
    
    %path taken by each channel
    path = zeros(K,Test_Limit);
    
    %distance history
    dist = zeros(K,Test_Limit);
    
    %indicies where good channels were found
    good_decision_trial = zeros(1,good);
    
    %running sum of ones
    d_k = zeros(1,K);
    
    %number of measurements each channel k has recieved
    m_k = ones(1,K);
    
    %The first minimum out threshold distance
    e_k = zeros(1,K);
    
    %the current cost
    cur_cost = zeros(1,K);
    
    %initalize cost and threshold distance
    for k = 1:1:K
        e_k(k) = min(abs(T_U(1) - d_k(k)),abs(T_L(1) - d_k(k)));
        cur_cost(k) = e_k(k);
    end
    
    %What decision was made for channel K
    decision = zeros(1,K);
    
    %when was it made (what global index)
    decided_at = zeros(1,K);
    
    %pick_sequence
    pick_sequence = zeros(1,Test_Limit);
        
    %write down in an error was made
    errors = zeros(1,K);
    
    for j = 1:1:Test_Limit
        %if All decision have been made stop this run
        if sum(decision > 0 ) == K
            break;
        end
        
        %update the cost
        for k = 1:1:K
            if decision(k) ~= 0
                %if a decision was made cost is infinity, this channel needs
                %no more measurements
                cur_cost(k) = max_cost;
            else
                %update the distance and cost
                U_dist = pdist([m_k(k), T_U(m_k(k)); m_k(k),d_k(k)],'euclidean');
                L_dist = pdist([m_k(k),T_L(m_k(k)); m_k(k),d_k(k)],'euclidean');
                e_k(k) = min(U_dist,L_dist);
                cur_cost(k) = e_k(k) + (lambda * m_k(k)) ;
            end
        end
        
        %count the deadline values
        if j == deadline
            deadline_comp_full(i) =  sum(decision > 0);
            sprintf('Number at deadline is %d',deadline_comp_full(i))
            
            %determining all the errors that happened upto the deadline
            dead_bin_p_act = bin_p_act(decision > 0);
            dead_decision = decision(decision > 0);
            dead_errors = zeros(1,length(dead_bin_p_act));
            for k = 1:1:length(dead_bin_p_act)
                if dead_bin_p_act(k) ~= dead_decision(k)
                    if dead_bin_p_act(k) == 2
                        %if the proper bin was 2, and we did not say 2, that is a type 2 error
                        dead_errors(k)=2;
                    else
                        %other wise it was a type 1 error (False rejection of null %hypothesis)
                        dead_errors(k)=1;
                    end
                end
            end
            dead_type1_rate = sum(dead_errors == 1) / length(dead_bin_p_act);
            dead_type2_rate = sum(dead_errors == 2) / length(dead_bin_p_act);
            Dead_Trial_error(i,:) = [dead_type1_rate,dead_type2_rate];
            sprintf('Dead trail error rates are alpha = %f, beta = %f ',Dead_Trial_error(i,1),Dead_Trial_error(i,2))
        end
        
        for k = 1:1:K
            %check for decisions
            if decision(k) == 0
                if d_k(k) > T_U(m_k(k))
                    %above the upper boundary is region 3
                    decision(k) = 3;
                    decided_at(k) = j;
                elseif d_k(k) < T_L(m_k(k))
                    %below the lower boundary is region 1, the "good"
                    %region since this means p is low
                    decision(k) = 1;
                    decided_at(k) = j;
                    good_decision_trial(find(good_decision_trial == 0, 1, 'first')) = j;
                elseif T_C_L(m_k(k)) < d_k(k) && d_k(k) < T_C_U(m_k(k))
                    %In between the center boundaries means we're in the
                    %second region
                    decision(k) = 2;
                    decided_at(k) = j;
                end
            end
        end
        
        %pick the channel smallest distance
        cur_pick = pick(cur_cost);
        pick_sequence(j) = cur_pick;
        %sample that channel
        d_k(cur_pick) = d_k(cur_pick) + random('bino',1,p_act(cur_pick));
        %increment it's meausrement count
        m_k(cur_pick) = m_k(cur_pick) + 1;
      
        cost(:,j) = cur_cost';
        path(:,j) = d_k';
        dist(:,j) = e_k';
    end
    
    %in case it finishes before the deadline
    if j < deadline
        deadline_comp_full(i) =  sum(decision > 0 );
        sprintf('Number at deadline is %d',deadline_comp_full(i))
    end
        
    %store decision times and compute errors
    for k = 1:1:K
        %flag the error types, the null hypothesis is that the paremter lives on
        %an edge (bin 1 or 3), if the bind don't match we have an error
        if bin_p_act(k) ~= decision(k)
            if bin_p_act(k) == 2
                %if the proper bin was 2, and we did not say 2, that is a type 2 error
                errors(k)=2;
            else
                %other wise it was a type 1 error (False rejection of null %hypothesis)
                errors(k)=1;
            end
        end
        
        if decided_at(k) == 0
            %this channel never completed, decison time is infinity
            full_decision_time(k,i) = Inf;
        else
            full_decision_time(k,i)= decided_at(k);
        end
    end
    
    type1_rate = sum(errors == 1) / K;
    type2_rate = sum(errors == 2) / K;
    Trial_error(i,:) = [type1_rate,type2_rate];
    
    %how long did it take to find 4 good ones.
    good_decision_time(i) = max(good_decision_trial);
    
end

comp_full = max(full_decision_time);
E_n = mean(comp_full);
E_n_good = mean(good_decision_time);
Dead_comp = mean(deadline_comp_full);
Error_rate = mean(Trial_error,1);
Dead_Error_rate = mean(Dead_Trial_error,1);

%store the results for comparison
save('./CC128N100L.mat');