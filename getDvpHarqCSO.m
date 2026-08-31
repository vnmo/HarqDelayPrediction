function dvp = getDvpHarqCSO(d, delta, pHarq, f, T, M, Qmax, CSO, zeta)
%GETDVPHARQCSO Proposed model: total DVP approximation.
% Analytical call path: getWaitPmfCSO -> getQLenPmfCSO for waiting and
% auxiliary overflow; getServiceDvpCSO for the remaining service budget.
% This function does not call the queueing simulation simulate_queue_IR_cso.
% d and T use the same time unit; d may be an array of deadlines.
% delta is the raw feedback delay in slots; zeta is the raw decoding delay.
% Each delay is rounded up separately by getServiceDvpCSO.
% The legacy eight-input call uses zeta=1, hence RTT=ceil(delta)+2.
% pHarq contains conditional PERs. CSO=1/xi, or zero without blockage.
% Drops, waiting beyond the deadline, and service violations are combined
% without counting the admission probability twice. Discards are violations.

if nargin < 9
    zeta = 1;
end
validateattributes(d, {'numeric'}, {'real','nonnan'});
validateattributes(T, {'numeric'}, {'scalar','real','finite','positive'});

% Validate the service parameters even when d is empty.
getServiceDvpCSO([], zeta, delta, pHarq, M, CSO);
[waitpdf, pDrop] = getWaitPmfCSO(pHarq, f, M, Qmax, CSO);
waitingSlots = 0:numel(waitpdf)-1;
dvp = zeros(size(d));
for index = 1:numel(d)
    residualSlots = d(index) / T - waitingSlots;
    exceedsDeadline = waitingSlots > floor(d(index) / T);
    waitingTail = sum(waitpdf(exceedsDeadline));
    serviceViolation = getServiceDvpCSO(residualSlots(~exceedsDeadline), ...
        zeta, delta, pHarq, M, CSO);
    admittedViolation = waitingTail ...
        + sum(waitpdf(~exceedsDeadline) .* serviceViolation);
    dvp(index) = pDrop + (1-pDrop) * admittedViolation;
end
end
