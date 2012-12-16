%A script to plot the V Mask of a sampleing given some parameters choices. 

%P = random('unif',0,1);
P = 0.85;
N = 10;
mu_0=.50;
Data = Sample(P,N);

SampSum = zeros(size(Data));
for i = 1:1:length(Data(1,:))
    SampSum(i)=sum(Data(1,1:1:i));
end

p_hat = (1/N) .* SampSum;

S_i = zeros(size(p_hat));

for j = 1:1:length(p_hat)
    S_i(j)=max(p_hat(j) - mu_0,0);
end

plot(1:1:N,S_i,1:1:N,p_hat);
axis([1 N 0 1])
legend('S_i','Estimate of p');
title(sprintf('H_0: p = %f; P actual = %f',mu_0,P));
