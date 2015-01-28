% Generates an 2 test simulated operator Curve These
% can be compared in a plot. Also Simulated  E(n), p_miss and p_fa
%vs p_range

%needed when running on linux.
%addpath '../lib'

%Multi Threading
matlabpool open;

%% Non SPRT params
Trials = 10000;

%Here we let the test limit be arbitrarly large as we want to figure out
%when the test actually completes. We use a break to short circut the run
%if it completes before the limit. This here purely to avoid infinite
%loops.
Test_Limit = 10000;

%range of P's that will be simulated.
p_range = 0:0.01:1;

%% SPRT params

%width of indiffrence region
delta = 0.125;

%what portion of the indiffrence is deidicated to beta
split = 0.5;

%Outer Bind width
p_prime = 0.2;

%bounds on alpha effective and beta effective
alpha = 0.01;
beta = 0.05;

%Upper and lower sequential test thresholds (aproximated)
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%% left side hypothesis threhold line computation
%tuneing parameters - center and width
p_prime_ls = p_prime;

%boundaries of accept and reject region
p_0_ls = p_prime_ls - (split * delta);
p_1_ls = p_prime_ls + ((1 - split) * delta);

%generate test lines
[L_1_ls,L_0_ls] = GenTestLines(alpha,beta,p_1_ls,p_0_ls,Test_Limit);

%% right side hypothesis threhold line computation
%tuneing parameters - center and width
p_prime_rs = 1 - p_prime_ls;

%boundaries of accept and reject region
p_0_rs = p_prime_rs - (split * delta);
p_1_rs = p_prime_rs + ((1 - split) * delta); 
%generate test lines
[L_1_rs,L_0_rs] = GenTestLines(alpha,beta,p_1_rs,p_0_rs,Test_Limit);

%% Recentering for fairness
Outer_lower_boundary = L_0_ls - 1 ;
Outer_upper_boundary = L_1_rs - 1;
Inner_lower_boundary = L_1_ls - 1;
Inner_upper_boundary = L_0_rs - 1 ;

%% Storage Vars

%Simulated operator curve P(H_0 | p) - Here H_0 is declare the parameter in
%an outer bin where outer bins are p > Outer_upper_boundary or
%p < Outer_lower_boundary. 

L = zeros(1,length(p_range));

%Mean number of samples required to complete said test
E_n = zeros(1,length(p_range));

%Probablity we declared H_1 when H_0 was true, declared bin 2, but really
%bin 1 or 3 was correct - Beta
p_miss = zeros(1,length(p_range));

%Probablity we declared H_0 when H_1 was true, declared bin 1 or 3, but bin
%2 was correct - Alpha
p_fa = zeros(1,length(p_range));

%% Begin Simulation - P loop 

for p_index = 1:1:length(p_range)

    sprintf('On p %f',p_range(p_index))
    %Decisions 1, 2, or 3 for the respective bin
    decision = zeros(1,Trials);
    
    %What was the index that the decision was made
    decided_at = zeros(1,Trials);
    
    %Was H_0 declared
    declared_H_0 = zeros(1,Trials);
    
    %Was there  miss error?
    miss = zeros(1,Trials);
    
    %Was there a false alarm error?
    fa = zeros(1,Trials);

    %% P computed Values
    
    % determine what Bin P belongs to
    if p_range(p_index) > p_prime_rs
        p_bin = 3;
    elseif p_range(p_index) < p_prime_ls
        p_bin = 1;
    else
        p_bin = 2;
    end
    
    %% Trial Loop
    parfor trial_index = 1:Trials
    %for trial_index = 1:1:Trials
   %     sprintf('On Trail %d',trial_index)
        %The running some of ones d_k
        running_sum = 0; 
        
        
        for Sample_index = 1:1:Test_Limit
            %% The actual SPRT Loop
            if running_sum > Outer_upper_boundary(Sample_index)
                decision(trial_index) = 3;
                break;
            elseif running_sum < Outer_lower_boundary(Sample_index)
                decision(trial_index) = 1;
                break;
            elseif Inner_lower_boundary(Sample_index) < running_sum && running_sum < Inner_upper_boundary(Sample_index)
                decision(trial_index) = 2;
                break;
            else
                running_sum = running_sum + random('bino',1,p_range(p_index));
            end 
        end
        
        %% Data collection / Error checking after the test
        %Length of the sequence
        decided_at(trial_index) = Sample_index;
        
        %Check for H_0 declaration
        if decision(trial_index) == 3 || decision(trial_index) == 1
            declared_H_0(trial_index) = 1;
        end
        
        %Check for Errors - If the bins don't agree
        if decision(trial_index) ~= p_bin
            if p_bin == 2
                % and if the actaul bin was 2, we made a false alarm error
                fa(trial_index) = 1;
            else
                % Other wise we made a miss error
                miss(trial_index) = 1;
            end
        end
    end
    
    %% Store The means for this value of the parameter    
    L(p_index) =  mean(declared_H_0);
    E_n(p_index) =  mean(decided_at);
    p_miss(p_index) =  mean(miss);
    p_fa(p_index) =  mean(fa);
end

matlabpool close;

%% Save data in file
%store the results for comparison
fname = strcat('OC2TestSim-',datestr(now, 'mmmdd.yyyy-HH.MM'),'.mat');
save(fname);
