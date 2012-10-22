function [ Index ] = Threshold(cusum, H, twiddle )
%Takes the tabular cusum and returns the first value above the threshold
%   cusum is a 2 X N Vector, the result of tabulur cusum on the RV
%   H is the threshold that We are looking to get past
%   twiddle is the small allownace to be below the threhold
Uindex = find(cusum(1,:) > (H - twiddle), 1, 'first');
Lindex = find(cusum(2,:) < (-1 * H + twiddle), 1, 'first');
if isempty(Uindex) && isempty(Lindex)
    Index = length(cusum);
elseif isempty(Uindex)
    Index = Lindex;
elseif isempty(Lindex)
    Index = Uindex;
else
    Index = min(Uindex,Lindex);
end
end

