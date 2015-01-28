function [ Samples ] = Sample1dim( P, N )
%SAMPLE Generates vectors of size 1 X N samples from a bernouli RV with
%reinveted the wheel again.

Samples = random('bino',1,P,[1,N]);

