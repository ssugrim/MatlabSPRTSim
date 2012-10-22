function [ L ] = leading( data )
%LEADING Returns the number of ones in the front of data an 1 x N vector.
% Returns the number of Leading ones in data. Computes
% the bitwise AND between a the subvector and a string of ones. If the
% sum matches the sub vector length, then the AND was posotive, and the sub
% vector was all ones. 
L = 0;
len = length(data);
while (L+1 <= len) && (sum(data(1:1:L+1) & ones(1,L+1)) == L+1)
    L = L+1;
end
end

