%Computes the number of required samples to get the sample variance below a
%certain limit. 

%number of channels
K = 5;

Test_Limit = 10000;

p_act = [0.1, 0.25, 0.5,  0.75, 0.9 ];

dat = zeros(K,Test_Limit);

for k = 1:1:K
    dat(k,:) = Sample1dim(p_act(k),Test_Limit);
end

run_expectation = zeros(K,Test_Limit);
run_variance= zeros(K,Test_Limit);

for n = 1:1:Test_Limit
    for k = 1:1:K
        run_expectation(k,n) = sum(dat(k,1:n))/ n;
        run_variance(k,n) = sqrt(sum((dat(k,1:n) - run_expectation(k,n)).^2)) / n;
    end
end

set(0,'DefaultAxesFontSize',14)
set(gca, 'FontSize', 14)

plotStyle = {'b--','g--','g-','r-','b-'};
figure(1)
hold on;

for k = 1:1:K
plot(1:1:Test_Limit,run_variance(k,:),plotStyle{ mod(k,length(plotStyle))+1},'LineWidth',2);
legendInfo{k} = sprintf('p_{actual} = %f', p_act(k));
end
legend(legendInfo);
%plot(1:1:Test_Limit,run_variance(2,:),':b');
axis([10 500 0 0.2]);
line([10,500],[0.04,.04],'color','k','LineWidth',2);
xlabel('Samples (j)');
ylabel('Sample Variance \frac{1}{j}\sum_{i=1}^{j}(x_i - m_j)');
title('Sample variance for two channels as a function of sample count');
%c1leg = sprintf('Channel 1, p_{actual} = %f', p_act(1));
%c2leg = sprintf('Channel 2, p_{actual} = %f', p_act(2));
%legend(c1leg,c2leg);