function [ outvector ] = count_ones(invector)
%COUNT_ONES Replaces each postion of invector with a running total of ones
%apparently I reinvented the wheel
outvector = cumsum(invector);
end

