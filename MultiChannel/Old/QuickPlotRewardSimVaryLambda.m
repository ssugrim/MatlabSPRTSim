%Computes Mean number of samples spent as a function of lambda

mean_samples = zeros(length(lambda_range),length(Channels));
Trunc_max = 100;

for k = 1:1:length(Channels)
    for i = 1:1:length(lambda_range)
        
        numerator = round(squeeze(histogram(k,i,1:Trunc_max))');
        denominator = sum(numerator);
        sprintf('k =%d, i=%d,denom = %d',k,i,denominator)
        mean_samples(i,k) = sum((numerator ./ denominator) .* (1:1:Trunc_max));
    end
end

mean_samples(isnan(mean_samples))=0;

plotStyle = {'b+-','bo-','b*-','bs-','bd-','b^-','bv-','b<-','b>-','r+-','ro-','r*-','rs-','rd-','r^-','rv-','r<-','r>-'};
hold;

legendInfo = {};

count = 1;
for k = [1,4,8,12,14]
    plot(lambda_range,mean_samples(:,k),plotStyle{mod(k,length(plotStyle))+1},'LineWidth',2);
    legendInfo{count} = sprintf('N = %d',Channels(k));
    count = count + 1;
end

legend(legendInfo);

xlabel('\lambda - Measurement Cost');
ylabel('Mean number of samples given to a channel while searching');
title({'Mean samples spent per channel to discover 1 good channel';'As \lambda increases this value converges,';'reguardless of the size of the pool'})

