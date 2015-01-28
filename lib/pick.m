function [ res ] = pick(vec)
%Randomly picks one of the lowest values of vector. Returns the index of
%that value.
first = find(vec == min(vec),1,'first');
res = first + floor(random('unif',0,length(find(vec == min(vec)))));
end

