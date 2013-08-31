function [ L_1, L_0 ] = GenTestLines( alpha, beta, p_1,p_0, length)
%GENTESTLINES Returns L_1 and L_0 vectors of the upper and lower boundary lines

%Given the OC parameters. Input Arts are:
%A, B - Sprt thresholds derived from \alpha and \beta
%p_zero, p_one - The fixed points of the OC curve, where L(p_zero) =
%\alpha and L(p_one)=beta
%length is the length of the vector we want as output

%independant sample indiex
m = 0:1:length-1;

%fundemental constants
A = (1 - beta) / alpha;
B = beta / ( 1 - alpha);

%convinent names for probablity ratios
p_upper = (1 - p_1) / (1 - p_0);
p_lower = (p_1) / (p_0);

%Graphical test intercepts and slope
h_0 = (log(B)) / (log(p_lower) - log(p_upper));
h_1 = (log(A)) / (log(p_lower) - log(p_upper));
s = (log(1 / p_upper)) / (log(p_lower) - log(p_upper));

%Graphical test Paramters and Boundaries.  Purely for grpahing examples, not used in simulation.
L_1 = (s .* m) + (h_1 .* ones(size(m)));
L_0 = (s .* m) + (h_0 .* ones(size(m)));
end

