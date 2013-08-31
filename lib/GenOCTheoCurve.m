function [ p,L,E_n ] = GenOCTheoCurve( alpha, beta,p_one, p_zero, gran )
%GENOCTHEOCURVE Generates a theoretical OC curve from parameters
%   Detailed explanation goes here
%dummy parameter
h = -500:gran:500;
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%Theoretical operator curve second coordinate (parametric) - L(p)
L = ((A.^h) - ones(size(h))) ./ (A.^h - B.^h);

%convinent names for probablity ratios
p_upper = (1 - p_one) / (1 - p_zero);
p_lower = (p_one) / (p_zero);

%Theoretical operator curve first coordinate (parametric) - p
p = (ones(size(h)) - p_upper.^h) ./ (p_lower.^h - p_upper.^h);

%Theoretical Expected sequence length (E_p(n))
E_n = ((L .* log(B)) + ((ones(size(h)) - L) .* log(A))) ./ ((p .* log(p_lower)) + ((ones(size(h)) - p) .* log(p_upper)));

end

