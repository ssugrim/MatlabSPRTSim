%An early script that found number of leading ones from the samples
%TODO: fix outdated function names
K = 20;
P = random('unif',0,1,[1,K]);
N = 10;
Data = Sample(P,N);
for i = 1:1:size(Data,1)
    leaders(i) = leading(Data(i,:));
end
Numof1 = transpose(sum(Data,2));
Measurements = length(find(leaders == 0)) + sum(leaders)