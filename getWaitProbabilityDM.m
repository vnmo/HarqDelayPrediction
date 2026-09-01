function [prob] = getWaitProbabilityDM(pdfPoint, queLength, pVector, M, tau)

    k = pdfPoint;   % total waiting time in slots
    q = queLength;

    tau = max(1, round(tau));

    if mod(k, tau) ~= 0
        prob = 0;
        return
    end

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