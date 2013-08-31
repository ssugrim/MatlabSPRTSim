function [ res ] = pick(vec)
%Pick's the lowest value of a Nx1 matrix, if equal pick randomly. 
%Returns index
if all(vec == min(vec))
    res = randi(length(vec),1,1);
else
    res = min(find(vec == min(vec)));
end

end

