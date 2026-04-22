function [prob] = getWaitProbabilityDM(pdfPoint, queLength, pVector, M, tau)
%GETWAITPROBABILITY
% Probability that q packets require total waiting time k slots,
% when each HARQ attempt consumes tau slots.
%
% Inputs:
%   pdfPoint  = k (total waiting time in slots)
%   queLength = q
%   pVector   = HARQ failure probabilities
%   M         = max number of attempts
%   tau       = slots consumed per attempt (positive integer)

    k = pdfPoint;   % total waiting time in slots
    q = queLength;

    % tau must be a positive integer
    tau = max(1, round(tau));

    % If total waiting time is not a multiple of tau,
    % it cannot be represented by integer number of attempts.
    if mod(k, tau) ~= 0
        prob = 0;
        return
    end

    % Convert slot-time into effective attempt-count time
    kEff = k / tau;

    if (kEff < q || kEff > M*q)
        prob = 0;
        return
    elseif length(pVector) == M
        prob = getprob_IR(kEff, q, pVector, M, M);
        return
    elseif length(pVector) == 1
        pVector = ones(1, M) * pVector;
        prob = getprob_IR(kEff, q, pVector, M, M);
        return
    else
        error("incompatible Probability Vector for IR");
    end
end