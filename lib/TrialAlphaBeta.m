function [ type_1, type_2] = TrialAlphaBeta( K, act_bin_k, bin_k)
%TRIALALPHABETA  % computes alpha and beta for a single trial

%% Nameing conventions
%K = # of channels

%bin_k = The bin a channel was placed in - Zero means the channel was not
%measured, other wise it should be valued 1,2,3 for the respective bin. 
% 1 is the bing for  \hat{p} < Lower Threshold, 
% and 3 is the bin for \hat{p} > Upper threshold

%act_bin_k = The actual bin the channel belongs to. None of these should be
%empty.

% the null hypothesis is the channel belongs in the outer bin
% If the correct bin is 2, but we said 1 or 3, this is a type 2 error (BETA). Fasle appceptance of the null hypothesis.
% IF the correct bin is 1 or 3  but we say 2, this is a type 1 error (ALPHA), False  rejection.


%% classify the errors made

errors = zeros(1,K);

%Compare Bin Values
for k = 1:1:K
    %Dissagreement between something I measured and the actual.
    if bin_k(k) ~= 0
        if act_bin_k(k) ~= bin_k(k)
            if act_bin_k(k) == 2
                %if the proper bin was 2, and we did not say 2, that is a type 2 error
                %we've falsed accepted the null hypothesis
                errors(k)=2;
            else
                %other wise it was a type 1 error (False rejection of null %hypothesis)
                errors(k)=1;
            end
        end
    end
end

%% compute alpha and beta
outer_bins = sum(bin_k == 1) + sum(bin_k == 3);
measured = sum(bin_k ~= 0);

%error rate for those that were measured
type_1 = sum(errors == 1) / max(measured - outer_bins,1);

type_2 = sum(errors == 2) / max(outer_bins,1);
end

