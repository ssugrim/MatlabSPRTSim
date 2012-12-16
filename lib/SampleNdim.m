function [ Samples ] = SampleNdim( P, N )
%SAMPLE Generates vectors of size length(P) X N samples from a bernouli RV with
% parameter given by P(i), the input vector. 
% Returns an Array of Sample sequneces of Length N.
% P is a 1 x K vector of parameters for the bernouli RV. 
% Where K is the number of sequences you want to generate.
Samples = zeros(size(P,2),N);
for i = 1:1:size(P,2)
     Samples(i,:) = random('bino',1,P(i),[1,N]);
 end
end

