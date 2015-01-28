% Generates an 2 test simulated operator Curve These
% can be compared in a plot. Also Simulated  E(n), p_miss and p_fa
%vs p_act

%needed when running on linux.
%addpath '../lib'

Trials = 1000;
K = 2;
Test_Limit = 5000;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;
delta = 0.05;

p_prime_1 = 0.3;
p_prime_2 = 1 - p_prime_1;

%threshold lines
[T_11u,T_10l] = GenTestLines(alpha,beta,p_prime_1,delta,Test_Limit);
[T_21u,T_20l] = GenTestLines(alpha,beta,p_prime_2,delta,Test_Limit);

decision = zeros(1,Trials);
decided_at = zeros(1,Trials);

p_act = 0:0.01:1;
L = zeros(1,length(p_act));
E_n = zeros(1,length(p_act));
p_miss = zeros(1,length(p_act));
p_fa = zeros(1,length(p_act));
for j = 1:1:length(p_act)
    p_act(j)
    fast_decision = zeros(1,Trials);
    fast_decided_at = zeros(1,Trials);
    decision = zeros(1,Trials);
    decided_at = zeros(1,Trials);
    declared_H_0 = zeros(1,Trials);
    
    miss = zeros(1,Trials);
    fa = zeros(1,Trials);
    
    for i = 1:1:Trials
        dat = Sample1dim(p_act(j),Test_Limit);
        runs = count_ones(dat);
        res_1 = SeqThreshold(T_11u,runs,T_10l);
        res_2 = SeqThreshold(T_21u,runs,T_20l);
        if res_1(1) < res_1(2) && res_2(1) < res_2(2)
            decision(i) = 3;
            declared_H_0(i) = 1;
            decided_at(i) = max(res_1(4),res_2(4));
            if p_prime_1 < p_act(j) && p_act(j) < p_prime_2
                [p_prime_1,p_act(j),p_prime_2]
                fa(i) = 1;
            end
        elseif res_1(3) > res_1(2) && res_2(3) > res_2(2)
            decision(i) = 1;
            declared_H_0(i) = 1;
            decided_at(i) = max(res_1(4),res_2(4));
            if p_prime_1 < p_act(j) && p_act(j) < p_prime_2
                [p_prime_1,p_act(j),p_prime_2]
                fa(i) = 1;
            end
        elseif isinf(res_1(4)) || isinf(res_2(4))
            error('Infinite Length test encounterd, Set a bigger limit');
        else
            decision(i) = 2;
            decided_at(i) = max(res_1(4),res_2(4));
            if p_act(j) < p_prime_1 ||  p_prime_2 < p_act(j)
                [p_prime_1,p_act(j),p_prime_2]
                miss(i) = 1;
            end
        end
    end
    L(j)= mean(declared_H_0);
    E_n(j) = mean(decided_at);
    p_miss(j) = mean(miss);
    p_fa(j) = mean(fa);
end

save('./OC2TestSim.mat');
