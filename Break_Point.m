%A script to plot the the ARL VS mu;

mu_0 = 0.5;
mu = 0:.05:1;
N = 100;
H = 3;

%setting these to zero generates the original C_i
twiddle = 0;
K = 0;
ustart = 0;
lstart = 0;

MAX_dev = zeros(size(mu));
trial = 50;

for i = 1:1:length(mu)
    trial_value = zeros(1,trial);
    for t = 1:1:trial
        samples = Sample1dim(mu(i),N);
        cusum = Gen_tab_x(samples,mu_0,K,ustart,lstart);
        %here we are only intrested in the achieved maximum
        trial_value(t) = max(max(cusum(1,:),max(cusum(2,:))));
    end
    MAX_dev(i) = mean(trial_value);
end

save('./Break_data.mat');