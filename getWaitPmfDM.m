function [probDistr] = getWaitPmfDM(pVector, f, M, Qmax, tau)

    tau = max(1, round(tau));

    qLenPmf = getQLenPmfDM(pVector, f, M, Qmax, tau);

    p_emptyQ = qLenPmf(1);
    qLenPmf(1) = [];

    Lq = length(qLenPmf);

    % max delay support in slots
    maxDelay = M * tau * Lq;

    % index 1 <-> delay 0
    waitDistribution = zeros(1, maxDelay + 1);

    for qLen = 1:Lq
        condWaitDistribution = zeros(1, maxDelay + 1);

        % minimum delay for qLen packets ahead
        for delay = qLen * tau : maxDelay
            wdensity = getWaitProbabilityDM(delay, qLen, pVector, M, tau);
            condWaitDistribution(delay + 1) = wdensity;
        end

        condWaitDistribution = condWaitDistribution / sum(condWaitDistribution + (sum(condWaitDistribution)==0));
        waitDistribution = waitDistribution + qLenPmf(qLen) * condWaitDistribution;
    end

    % overwrite delay 0 with empty-queue mass
    probDistr = waitDistribution;
    probDistr(1) = p_emptyQ;
end