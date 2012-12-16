function [res ] = SeqThreshold( uvec,dvec,lvec )
%SEQTHRESHOLD Returns the [Upper Value, Data value, Lower vale, Index] of the threshold crossing
%   Requires 3 vectors, uvec, dvec, and lvec each of size Nx1. Uvec is the 
%upper limit vecotor Line (L_1) and lvec is the lower anagloue. dvec is the
%data vector. Will loop through indicies until it finds the corssing, and
%return the upper, Lower and data values, and index that occured. 

lcomp = dvec > lvec;
ucomp = uvec > dvec;

lind = find(lcomp==0, 1, 'first');
uind = find(ucomp==0, 1, 'first');

ind = min([lind,uind]);
res = [uvec(ind),  dvec(ind), lvec(ind), ind];


end

