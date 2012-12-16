function [ outvector ] = count_ones(invector)
%COUNT_ONES Replaces each postion of invector with a running total of ones
outvector = zeros(size(invector));
for i = 1:1:length(invector)
    outvector(i)=sum(invector(1:i));
end

end

