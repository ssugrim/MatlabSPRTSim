% Generates data and plots the CUSUM under given parameters

mu = 0.89;
mu_0 = 0.95;
%sample length
N = 10000;
%slack value
K=.2;
%head start
ustart=0;
lstart=0;
%graph max
Y_max = 5;

%threshold to graph
H = 3;

%generate the samples and CUSUM via my functions
samples = Sample1dim(mu,N);
cusum = Gen_tab_x(samples,mu_0,K,ustart,lstart);

%plot and hold for other data
plot(cusum');
hold on;

%threhsold lines
topline = H * ones(1,N);
botline = -H * ones(1,N);
plot(0:N-1,topline,'red',0:N-1,botline,'black');

%plot formatting and labeling
legend('C^{+}_{i}','C^{+}_{i}','y=3','y=-3')
axis([0 3000 -Y_max Y_max])
xlabel('step index, i');
ylabel('C_i');
mu_0_str = sprintf('%5.2f',mu_0);
mu_str = sprintf('%5.2f',mu);
n_str = sprintf(' N = %6.0f,',N);
k_str = sprintf(' K = %5.2f',K);
title(strcat('Tabular CUSUM for \mu_0 =',mu_0_str,' and \mu =',mu_str,'; with parameters:',n_str,k_str));
hold off;