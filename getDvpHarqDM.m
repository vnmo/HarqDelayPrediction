function [dvp] = getDvpHarqDM(d, delta, pHarq, f, T, M, Qmax, tau)

    tau = max(1, round(tau));

    waitpdf = getWaitPmfDM(pHarq(:), f, M, Qmax, tau);

    dvp = 0;

    for delay = 0:length(waitpdf)-1
        k_d = min(M, floor(((d / T) - delay + delta) / (delta + tau)));

        if k_d <= 0
            failProb = 1;
        else
            % if pHarq is cumulative BLER, use failProb = pHarq(k_d)
            % if pHarq is conditional failure probability, use product
            failProb = prod(pHarq(1:k_d));
        end

        dvp = dvp + waitpdf(delay + 1) * failProb;
    end
end