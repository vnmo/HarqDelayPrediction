function [dvp] = getDvpHarq (d, delta, pHarq, f, T, M, Qmax)



waitpdf = getWaitPmf(pHarq(:), f, M, Qmax);
dvp = 0;
for k = 0:length(waitpdf) - 1
    k_d = min(M, floor(((d / T) - k + delta) / (delta + 2)));
    
    dvp = dvp + waitpdf(k + 1) * prod(pHarq(1:k_d));
end


end
