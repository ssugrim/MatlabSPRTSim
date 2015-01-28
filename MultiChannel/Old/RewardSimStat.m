%Generates "Data" for a 2 channel reward model. 

%for use during cmd line run in linux
%addpath '../lib'

%global Parameters

%number of channels (specifiying specific channel layout)
p_act = [0.0193,0.2113,0.0368,0.4656,0.2159,0.2251,0.1312,0.2975,0.1609,0.4347];
K = length(p_act);

%number of good channels we want to dig out
good = 4;

%dealine to count Completed (we'll let the test run to completeion, we
%merely count here).
deadline = 2000;

%penalitly coefficent;
lambda = 0.1;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

delta = 0.1;
p_prime = 0.3;
split = 0.3;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%The hard limt All tests must stop at (Lenght of our vectors, should be
%longer than eta, for some allowance).
Test_Limit = 10000;

%number of Trials
Trials = 10000;

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

full_decision_time = zeros(K,Trials);
part_decision_time = zeros(K,Trials);
good_decision_time = zeros(1,Trials);
good_available = zeros(1,Trials);
deadline_comp_full = zeros(1,Trials);
deadline_comp_part = zeros(1,Trials);

%Begin Simulation

for i = 1:1:Trials
    sprintf('On trail %d',i)
    %per trial channel Parameters (rand draws from uniform 0,1)
    
    % for Random Trials (commented to tests a specif vector)
   % p_act = rand(1,K);
    good_available(i) = sum(p_act < p_prime_ls);
    %Scenario Comparison for now.
    %p_act = [0.9, 0.6];
    
    runs = zeros(K,Test_Limit);
    dists = zeros(K,Test_Limit);
            
    %generate the samples for each channel
    for k = 1:1:K
        dat = Sample1dim(p_act(k),Test_Limit);
        runs(k,:) = count_ones(dat);
        dists(k,:) = min(abs(L_0_ls - runs(k,:)),abs(L_1_rs - runs(k,:)));
    end
    
    %compute the largest possible cost
    maxcost = max(max(dists)) + (lambda * Test_Limit);
    
    %Per trial Run Vars (can be used for plotting).
    cost = zeros(K,Test_Limit);
    clocks = zeros(K,Test_Limit);
    
    cur_clocks = ones(K,1);
    fast_decision = zeros(K,1);
    fast_decided_at = zeros(K,1);
    decision = zeros(K,1);
    decided_at = zeros(K,1);
        
    cur_pick = pick(dists(:,1));
    cur_clocks(cur_pick)= cur_clocks(cur_pick) + 1;
    cur_cost = dists(:,1);
    cost(:,1)= cur_cost;
    clocks(:,1) = cur_clocks;
    good_decision_trial = zeros(1,good);
    
    for j = 2:1:Test_Limit
        
        %count the deadline values 
        if j == deadline
          deadline_comp_full(i) =  sum(decision > 0 );
          deadline_comp_part(i) = sum(fast_decision > 0);
        end
        
        for k = 1:1:K
            %Check for Decision 
            if fast_decision(k) == 0
                %first make the fast decison, regions are numbered left to
                %right 1,2,3. Region 4 = (2,3) and region 5 (1,2)
                if runs(k,cur_clocks(k)) > L_1_ls(cur_clocks(k))
                    fast_decision(k) = 4;
                    fast_decided_at(k) = j;
                elseif runs(k,cur_clocks(k)) < L_0_rs(cur_clocks(k))
                    fast_decision(k) = 5;
                    fast_decided_at(k) = j;
                end
            elseif fast_decision(k) == 4 && decision(k) == 0
                %decide between (1,2)
                if runs(k,cur_clocks(k)) < L_0_rs(cur_clocks(k))
                    decision(k) = 2;
                    decided_at(k) = j;
                elseif runs(k,cur_clocks(k)) > L_1_rs(cur_clocks(k))
                    decision(k) = 3;
                    decided_at(k) = j;
                end
            elseif fast_decision(k) == 5 && decision(k) == 0
                %decide between (2,3)
                if runs(k,cur_clocks(k)) < L_0_ls(cur_clocks(k))
                    decision(k) = 1;
                    decided_at(k) = j;
                    good_decision_trial(find(good_decision_trial == 0, 1, 'first')) = j;
                elseif runs(k,cur_clocks(k)) > L_1_ls(cur_clocks(k))
                    decision(k) = 2;
                    decided_at(k) = j;
                end
            end
        end

        %pick the smallest distance
        cur_pick = pick(cur_cost);
        %update the clock
        cur_clocks(cur_pick)= cur_clocks(cur_pick) + 1;
        
        %update the reward
        for k = 1:1:K
            if decision(k) ~= 0
                %if a decision was made, zero reward
                cur_cost(k) = maxcost;
            else
                % current reward
                cur_cost(k) = dists(k,cur_clocks(k)) + (lambda * cur_clocks(k)) ;
            end            
        end
        %record the reward as a function of time
        cost(:,j)= cur_cost;
        clocks(:,j) = cur_clocks;
        
        if sum(decision > 0 ) == K
            break;
        end
    end
    
    %store decision times
    
    %full decision
    for k = 1:1:K
        if decided_at(k) == 0
            full_decision_time(k,i) = Inf;
        else
            full_decision_time(k,i)= decided_at(k);
        end
    end
    
    %partial decision
    for k = 1:1:K
        if fast_decided_at(k) == 0
            part_decision_time(k,i) = Inf;
        else
            part_decision_time(k,i)= fast_decided_at(k);
        end
    end
    
    good_decision_time(i) = max(good_decision_trial);
    
end

comp_full = max(full_decision_time);
comp_part = max(part_decision_time);
E_n = mean(comp_full);
E_n_good = mean(good_decision_time);


%store the results for comparison
save('./rewardsimdata.mat');