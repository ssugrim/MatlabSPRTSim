function [ Samples ] = Sample1dim( P, N )
%SAMPLE Generates vectors of size 1 X N samples from a bernouli RV with
% parameter given by P. 
% Returns an Array of Sample sequneces of Length N.
% P is a probablity between 0 and 1. 
Samples = zeros(1,N);
for i = 1:1:length(Samples)
     Samples(1,i) = random('bino',1,P);
 end
end

