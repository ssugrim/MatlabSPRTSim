%Generates Plots for a one run of a K channel reward model. 

%for use during cmd line run in linux
%addpath '../lib'

%global Parameters

%number of channels
K = 5;

%penalitly coefficent;
lambda = 0.10;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;
p_prime = 0.3;
delta = 0.1;
split = 0.3;


%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%The hard limt All tests must stop at (Lenght of our vectors, should be
%longer than eta, for some allowance).
Test_Limit = 10000;

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

%per trial channel Parameters (rand draws from uniform 0,1)
%p_act = rand(1,K)
p_act = [0.1, 0.25, 0.5, 0.75, 0.9];
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

%storage Vars
cur_clocks = ones(K,1);
fast_decision = zeros(K,1);
fast_decided_at = zeros(K,1);
decision = zeros(K,1);
decided_at = zeros(K,1);

%startup sequence
cur_pick = pick(dists(:,1));
cur_clocks(cur_pick)= cur_clocks(cur_pick) + 1;
cur_cost = dists(:,1);

for j = 2:1:Test_Limit
    %per run transient vars
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
end

% measurements  = zeros(K,Test_Limit);
% for k = 1:1:K
%     measurements(k,:) = count_ones(measured(k,:));
% end


plotStyle = {'r+-','ro-','r*-','rs-','rd-'};

fast_decision 
fast_decided_at
decision 
decided_at 

set(0,'DefaultAxesFontSize',20)
set(gca, 'FontSize', 20)

%The actual tests 
% figure(1)
% hold on;
% plot(L_1_ls,'color','k');
% legendInfo1{1} = sprintf('L_{1,1}');
% plot(L_0_ls,'color','g');
% legendInfo1{2} = sprintf('L_{1,0}');
% plot(L_1_rs,'color','b');
% legendInfo1{3} = sprintf('L_{2,1}');
% plot(L_0_rs,'color','m');
% legendInfo1{4} = sprintf('L_{2,0}');
% 
% for k = 1:1:K
%     plot(runs(k,:),plotStyle{ mod(k,length(plotStyle))+1});
%     legendInfo1{k+4} = sprintf('Channel %d, p_{act} = %f',k,p_act(k));
% end
% x_end = round(max(fast_decided_at))+10;
% y_max = round(max(runs(:,x_end)))+10;
% axis([1,x_end,0,y_max]);
% legend(legendInfo1)
% titlestr = sprintf('Sample Threshold,p` = %f, delta = %f, split = %f ',p_prime,delta,split);
% title(titlestr);
% xlabel('Time Index (j)');
% ylabel('Cost');

hold off;


%Cost per channel
plotStyle = {'b-','g-','k-','m-','r-'};
figure(1)
hold on;

for k = 1:1:K
    plot(cost(k,:),plotStyle{ mod(k,length(plotStyle))+1},'LineWidth',2);
    legendInfo2{k} = sprintf('Cost Channel%d, p_{act} = %f, decided at %f',k,p_act(k),decided_at(k));
end
max_at = round(max(decided_at));
max_rnd = max_at - mod(max_at, 100) + 50
axis([0,max_rnd,0,round( max(max(cost(cost < maxcost))))+10]);

legend(legendInfo2);
titlestr = sprintf('Cost per channel, lambda = %f ',lambda);
title(titlestr);
xlabel('Time Index (m)');
ylabel('Cost');
hold off;

%measurements 
figure(2)
hold on;

for k = 1:1:K
    plot(clocks(k,:),plotStyle{ mod(k,length(plotStyle))+1},'LineWidth',2);
    legendInfo3{k} = sprintf('Measurements Channel %d, p_{act} = %f, decided at %f',k,p_act(k),decided_at(k));
end
axis([0,max_rnd,0,round(max(max(clocks(:,1:max_rnd)))+10)]);
legend(legendInfo3);
titlestr = sprintf('Measurements per channel, lambda = %f ',lambda);
title(titlestr);

%line([max_at,max_at],[0,round(max(max(clocks(:,1:max_rnd)))+10)],'color','k');

xlabel('Time Index (m)');
ylabel('Measurements');
hold off;

%store the results for comparison
save('./rewardsimdata.mat');