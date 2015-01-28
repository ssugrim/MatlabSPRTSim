%Computes the number of required samples to get the sample variance below a
%certain limit.

%for use during cmd line run in linux
%addpath '../lib'

Test_Limit = 10000;
Trials = 100;

var_threshold = 0.01;

p_act = 0.01:0.01:0.99;
samples = zeros(size(p_act));

for k = 1:1:length(p_act)
    %per Trial Storage vars
    p_act(k)
    samp_trial = zeros(1,Trials);
    
    for j = 1:1:Trials
        %generate samples
        dat = Sample1dim(p_act(k),Test_Limit);
        
        %tmp storage vars
        run_expectation = zeros(size(dat));
        run_variance= zeros(size(dat));
        
        %compute variances
        for n = 1:1:Test_Limit
            run_expectation(n) = sum(dat(1:n))/ n;
            run_variance(n) = sqrt(sum((dat(1:n) - run_expectation(1:n)).^2)) / n;
        end
        
        %find # of samples
        max_comp = (run_variance == max(run_variance));
        max_ind =  find(max_comp==1, 1, 'first');
        run_var_comp = (run_variance(max_ind:end) < var_threshold);
        ind = find(run_var_comp==1, 1, 'first');
        if isempty(ind)
            samp_trial(j) = Inf;
        else
            samp_trial(j) = ind;
        end
    end
    
    %store mean
    samples(k) = mean(samp_trial);
end

save('./SimpleSample.mat');

% figure(1);
% hold on;
% plot(p_act, samples);
% xlabel('Paramter p_{actual}');
% ylabel('Expected numer of samples to complete test');
% line([0,1],[100,100],'color','r');
% line([0,1],[50,50],'color','r');
% titlestr = sprintf('Expected number of samples for variance <  %f', var_threshold);
% title(titlestr);
% 
