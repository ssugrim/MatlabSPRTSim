% Generates an 2 test simulated operator Curve These
% can be compared in a plot. Also Simulated  E(n), p_miss and p_fa
%vs p_range

%needed when running on linux.
%addpath '../lib'


%Here we let the test limit be arbitrarly large as we want to figure out
%when the test actually completes. We use a break to short circut the run
%if it completes before the limit. This here purely to avoid infinite
%loops.
Test_Limit = 35;

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

s = [0,0,1,0,0,0];
ideal_bot = zeros(1,length(L_0_ls) - length(s)+ 1);
ideal_top = ones(1,length(L_0_ls) - length(s) + 1);

f1 = figure(1);
hold on;
set(findall(f1,'type','text'),'fontSize',36,'fontWeight','bold','LineWidth',4);
plot(0:length(Outer_lower_boundary) - 1, Outer_lower_boundary,'g--*','LineWidth',3);
plot(0:length(Inner_lower_boundary) - 1, Inner_lower_boundary,'g--*','LineWidth',3);
plot(0:length(Outer_upper_boundary) - 1, Outer_upper_boundary,'b-s','LineWidth',3);
plot(0:length(Inner_upper_boundary) - 1, Inner_upper_boundary,'b-s','LineWidth',3);
plot((length(s) - 1) + (0:(length(ideal_bot)-1)),max(cumsum(s)) + cumsum(ideal_bot),'m-o','LineWidth',4);
plot((length(s) - 1) + (0:(length(ideal_top)-1)),(max(cumsum(s)) + cumsum(ideal_top) - 1) ,'c-o','LineWidth',4);
plot(0:length(s) - 1,cumsum(s),'r-o','LineWidth',4);
xlabel('Sample Index (m_n)');

