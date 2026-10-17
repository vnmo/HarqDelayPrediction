function dvp = getDvpHarqCSO(d, delta, pHarq, f, T, M, Qmax, CSO, zeta)
% Inputs: deadline d, feedback delay delta, conditional PER pHarq, arrival probability f, slot duration T, attempt limit M, maximum aggregate packet count Qmax, CSO, and decoding delay zeta.
% Output: DVP upper bound for the phase-averaged aggregate model.
% Theorem 1: overflow + admission*(waiting tail + finite-wait service term).
% This is not an analytical bound on a finite-process physical HARQ system.
% The waiting-delay function is exact within the stated aggregate model;
% the upper-bound inequality comes from the service-delay bound.
% The final min(1,...) retains the universal probability cap.

if nargin < 9
    zeta = 1;
end
validateattributes(d, {'numeric'}, {'real','nonnan'});
validateattributes(T, {'numeric'}, {'scalar','real','finite','positive'});

getServiceDvpCSO([], zeta, delta, pHarq, M, CSO);
[waitFunction, pDrop] = getWaitFunctionCSO(pHarq, f, M, Qmax, CSO);
waitingSlots = 0:numel(waitFunction)-1;
dvp = zeros(size(d));
for index = 1:numel(d)
    residualSlots = d(index) / T - waitingSlots;
    exceedsDeadline = waitingSlots > floor(d(index) / T);
    waitingTail = sum(waitFunction(exceedsDeadline));
    serviceViolation = getServiceDvpCSO(residualSlots(~exceedsDeadline), ...
        zeta, delta, pHarq, M, CSO);
    admittedViolation = waitingTail ...
        + sum(waitFunction(~exceedsDeadline) .* serviceViolation);
    dvp(index) = min(1, pDrop + (1-pDrop) * admittedViolation);
end
end
