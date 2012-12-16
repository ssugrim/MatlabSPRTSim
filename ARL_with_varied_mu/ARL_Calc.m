%A script to Compute the Values of the ARL for first crossing of Tabular
%Cusum with given parameters across various values of \mu. 

mu_0 = 0.5;
mu = 0:.05:1;
N = 10000;
H = 3;
twiddle = 0;
K = .2;
ustart = 1.5;
lstart = 1.5;
ARL = zeros(size(mu));
trial = 50;

for i = 1:1:length(mu)
    trial_value = zeros(1,trial);
    for t = 1:1:trial
        samples = Sample1dim(mu(i),N);
        cusum = Gen_tab_x(samples,mu_0,K,ustart,lstart);
        trial_value(t) = Threshold(cusum,H,twiddle);
    end
    ARL(i) = mean(trial_value);
end

save('./ARL_data.mat','ARL','mu_0','trial','N','H','K','ustart','lstart','mu');