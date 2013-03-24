function [ data ] = GenOCthreshold( P, s, h_0, h_1, deadline)
%GENOCPOINT Computes the crossing values for a given value of P, and the thereshold paramters
% P is the actual probablity of getting a one. s is the common slope of the
% the two thresholds. h_0 is the lower threhold intercept, and h_1 is the
% upper threhsold itnercept. deadline is the hard limit we want to stop at,
% we return all zeros if we nevre reach the deadline.

%data variable initilization

%number of ones
d = 0;

%ending index
ind = 0;

%upper limit at crossing
ul = 0;

%lower limit at corssing
ll = 0;

for i = 1:1:deadline
    %generate a sample from the distribution
    d = d + random('bino',1,P);
    
    %check ll
    if d < ((s * i) + h_0)
        ll = ((s * i) + h_0);
        ul = ((s * i) + h_1);
        ind = i;
        break;
    end
    
    %check ul
    if d > ((s * i) + h_1)
        ll = ((s * i) + h_0);
        ul = ((s * i) + h_1);
        ind = i;
        break;
    end
end

%return value
data = [ul, d, ll, ind];
end

