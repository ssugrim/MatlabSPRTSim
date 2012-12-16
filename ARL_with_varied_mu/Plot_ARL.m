%Generates the PLOT of ARL from ARL_data.mat generated ARL_Calc.m

load('I:\My Documents\Academics\research\Channel Surfing\Simulation\SampleGen\ARL_data.mat')
plot(mu,ARL);
axis([0 1 0 1000])
legend('ARL to detect shift');
xlabel('\mu actual');
str = sprintf('ARL over %5.2f trials',trial);
ylabel(str);
mu_str = sprintf('%5.2f',mu_0);
n_str = sprintf('%5.2f',N);
h_str = sprintf(', H = %5.2f',H);
k_str = sprintf(', K = %5.2f',K);
title(strcat('ARL to Detect \mu diffrent from H_0: \mu_0 = ',mu_str,'; with parameters: N = ',n_str,h_str,k_str));