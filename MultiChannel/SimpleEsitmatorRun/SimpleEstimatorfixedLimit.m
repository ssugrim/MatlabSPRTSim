%Computes the number of channels correcly categorised by a simple estimator
% that divides the number of samples allowed evenly among all the channels.

%number of channels
K = 256;

%Test Limit - Total number of samples I can allocate
Test_Limit = 1000;

T_L = 0.2;
T_R = 1 - T_L;

%sample allotment
sample_allot = floor(Test_Limit / K);

%p_act = [0.1, 0.25, 0.5,  0.75, 0.9 ];

%Number of Trials to run
Trials = 100;

Trial_error = zeros(Trials,2);

for t = 1:1:Trials 
    %Draw 64 Random channels
    p_act = rand(1,K);
    
    sprintf('On trail %d',t)
    
    %storage vars
    dat = zeros(K,sample_allot);
    hat_p = zeros(1,K);
    bin_p_hat = zeros(1,K);
    bin_p = zeros(1,K);
    errors = zeros(1,K);

    for k = 1:1:K
        %Sample the kth channel the allotted times
        dat(k,:) = Sample1dim(p_act(k),sample_allot);
        
        %compute the estimate
        hat_p(k) = sum(dat(k,:))/ sample_allot;
        
        %bin the estimates
        if hat_p(k) > T_R
            bin_p_hat(k) = 3;
        elseif hat_p(k) < T_L
            bin_p_hat(k) = 1;
        else
            bin_p_hat(k) = 2;
        end
        
        %bin the actual
        if p_act(k) > T_R
            bin_p(k) = 3;
        elseif p_act(k) < T_L
            bin_p(k) = 1;
        else
            bin_p(k) = 2;
        end
        
        %flag the error types, the null hypothesis is that the paremter lives on
        %an edge (bin 1 or 3), if the bind don't match we have an error
        if bin_p(k) ~= bin_p_hat(k)
            if bin_p(k) == 2
                %if the proper bin was 2, and we did not say 2, that is a type 2 error
                errors(k)=2;
            else
                %other wise it was a type 1 error (False rejection of null
                %hypothesis)
                errors(k)=1;
            end
        end
    end
    
    %compute the error rates
    
    type1_rate = sum(errors == 1) / K;
    type2_rate = sum(errors == 2) / K;
    
    %store this trials error rate
    Trial_error(t,:) = [type1_rate,type2_rate];
end

Dead_error_rate = mean(Trial_error);

fname = sprintf('SE%dT%d.mat',K,Trials);
save(fname);
% set(0,'DefaultAxesFontSize',14)
% set(gca, 'FontSize', 14)
% 
% plotStyle = {'b--','g--','g-','r-','b-'};
% figure(1)
% hold on;
% 
% for k = 1:1:K
% plot(1:1:Test_Limit,run_variance(k,:),plotStyle{ mod(k,length(plotStyle))+1},'LineWidth',2);
% legendInfo{k} = sprintf('p_{actual} = %f', p_act(k));
% end
% legend(legendInfo);
% %plot(1:1:Test_Limit,run_variance(2,:),':b');
% axis([10 500 0 0.2]);
% line([10,500],[0.04,.04],'color','k','LineWidth',2);
% xlabel('Samples (j)');
% ylabel('Sample Variance \frac{1}{j}\sum_{i=1}^{j}(x_i - m_j)');
% title('Sample variance for two channels as a function of sample count');
% %c1leg = sprintf('Channel 1, p_{actual} = %f', p_act(1));
% %c2leg = sprintf('Channel 2, p_{actual} = %f', p_act(2));
%legend(c1leg,c2leg);