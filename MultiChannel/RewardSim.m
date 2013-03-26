%Generates "Data" for a 2 channel reward model. 

%for use during cmd line run in linux
%addpath '../lib'

%global Parameters

%number of channels
K = 2;

%rweard coefficent
rho = 2;

%penalitly coefficent;
lambda = 1;

%Per channel limit eta (currently the point where the inner bounds of the
%tests intersect (deadline for an fast decision).
eta = 40;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%The hard limt All tests must stop at (Lenght of our vectors, should be
%longer than eta, for some allowance).
Test_Limit = 100;

%number of Trials
Trials = 1;

%%%%%%%%%% right side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_rs = 0.7;
delta_rs = 0.05;

%boundaries of accept and reject region
p_1_rs = p_prime_rs + delta_rs;
p_0_rs = p_prime_rs - delta_rs;

%generate test lines
[L_1_rs,L_0_rs] = GenTestLines(A,B,p_1_rs,p_0_rs,Test_Limit);

%%%%%%%%%% left side hypothesis %%%%%%%%%%%%%%%

%tuneing parameters - center and width
p_prime_ls = 0.3;
delta_ls = 0.05;

%boundaries of accept and reject region
p_1_ls = p_prime_ls + delta_ls;
p_0_ls = p_prime_ls - delta_ls;

%generate test lines
[L_1_ls,L_0_ls] = GenTestLines(A,B,p_1_ls,p_0_ls,Test_Limit);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Begin Simulation

for i = 1:1:Trials
    %per trial channel Parameters (rand draws from uniform 0,1)
   % p_act = rand(1,K);
    p_act = [0.9, 0.6];
    runs = zeros(K,Test_Limit);
            
    %generate the samples for each channel
    for k = 1:1:K
        dat = Sample1dim(p_act(k),Test_Limit);
        runs(k,:) = count_ones(dat);
    end
    
    %Per trial Run Vars (can be used for plotting).
    reward = zeros(K,Test_Limit);
    measured  = zeros(K,Test_Limit);
    
    fast_decision = zeros(K,1);
    fast_decided_at = zeros(K,1);
    decision = zeros(K,1);
    decided_at = zeros(K,1);
    
    for j = 1:1:Test_Limit
        %per run transient vars
        for k = 1:1:K
            %reward primitives
            dist_max = L_1_rs(Test_Limit) - L_0_rs(Test_Limit);
            dist = calc_dist(runs(k,j),L_1_rs(j),L_0_ls(j));
            spent = sum(measured(k,:)) / eta;
            reward(k,j) = (rho *(1 - (dist / dist_max))) - (lambda * spent);
                      
            if fast_decision(k) == 0
                %first make the fast decison, regions are numbered left to
                %right 1,2,3. Region 4 = (2,3) and region 5 (1,2)
                if runs(k,j) > L_1_ls(j)
                    fast_decision(k) = 4;
                    fast_decided_at(k) = j;
                elseif runs(k,j) < L_0_rs(j)
                    fast_decision(k) = 5;
                    fast_decided_at(k) = j;
                end
            elseif fast_decision(k) == 4 && decision(k) == 0
                %decide between (1,2)
                if runs(k,j) < L_0_rs(j)
                    decision(k) = 2;
                    decided_at(k) = j;
                elseif runs(k,j) > L_1_rs(j)
                    decision(k) = 3;
                    decided_at(k) = j;
                end
            elseif fast_decision(k) == 5 && decision(k) == 0
                %decide between (2,3)
                if runs(k,j) < L_0_ls(j)
                    decision(k) = 1;
                    decided_at(k) = j;
                elseif runs(k,j) > L_1_ls(j)
                    decision(k) = 2;
                    decided_at(k) = j;
                end
            end
            if decision(k) ~= 0
                %if no decision was made zero reward
                reward(k,j) = 0;
            end
        end
        if decision(k) == 0
            picked = find((reward(:,j)==max(max(reward(:,j)))));
            measured(picked,j) = 1;
        end
    end
end
measurements  = zeros(K,Test_Limit);
for k = 1:1:K
    measurements(k,:) = count_ones(measured(k,:));
end

%fast_decision 
%fast_decided_at
%decision 
decided_at 

figure(1)
hold on;
plot(L_1_ls,'color','r');
plot(L_0_ls,'color','g');
plot(L_1_rs,'color','b');
plot(L_0_rs,'color','m');
plot(runs(1,:),'+-','color','k');
plot(runs(2,:),'.-','color','k');
legend('Left side upper','Left side lower','Right side Upper','Right Side lower','Channel 1, p_{act} =0.9','Channel 2, p_{act} =0.6');
xlabel('Time Index (j)');
ylabel('Thresholds and Data');
y_lims = ylim;
y_lims(1) = 0;
ylim(y_lims);
hold off;

figure(2)
hold on;
plot(reward(1,:),'color','b');
plot(reward(2,:),'color','r');
legend('Reward Channel 1 p_{act} =0.9','Reward Channel 2 p_{act} =0.6')
xlabel('Time Index (j)');
ylabel('Reward');
hold off;

figure(3)
hold on;
plot(measurements(1,:),'color','g');
plot(measurements(2,:),'color','k');
legend('Measurements Channel 1 p_{act} =0.9','Measurements Channel 2 p_{act} =0.6')
xlabel('Time Index (j)');
ylabel('Measurements');
hold off;

%store the results for comparison
save('./rewardsimdata.mat');