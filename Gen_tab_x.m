function [ result ] = Gen_tab_x( sample, mean, K, ustart, lstart)
%GEN_TAB_X Generates a table directly from sample - mean
%   Returns a N x 2 matrix (where N = lenght of imput sample
%   vector). Calculates C^+_i and C^-_i from pg 404 of montogmery.
%   mean - hypothesised mean
%   sample - Nx1 vector of samples
%   K - Shift in mean we are willing to "tolerate"
%   ustart - Upper head start
%   lstart - Lower head start
result = zeros(2,length(sample));
for i = 1:1:length(sample)
    if (i-1)==0
        result(1,i) = max(0,(sample(1,i) - mean - K + ustart));
        result(2,i) = min(0,(sample(1,i) - mean + K + lstart));
        %result(1,i) = (sample(1,i) - mean - K + ustart);
        %result(2,i) = (sample(1,i) - mean + K + lstart);
    else
        result(1,i) = max(0,(sample(1,i) - mean - K + result(1,i-1)));
        result(2,i) = min(0,(sample(1,i) - mean + K + result(2,i-1)));
        %result(1,i) = (sample(1,i) - mean - K + result(1,i-1));
        %result(2,i) = (sample(1,i) - mean + K + result(2,i-1));
    end
 end
return
