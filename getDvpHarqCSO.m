function [dvp] = getDvpHarqCSO(d , delta, pHarq, f, T, M, Qmax,CSO)
%GETDVPHARQ Function to get the DVP with IR
%   Function description

% k_d = min(M, 1 + floor((d/T - 1) / (delta + 1)));
% dvp_IR_serv = prod(pHarq(1:k_d));
waitpdf = getWaitPmfCSO(pHarq(:), f, M, Qmax, CSO);
dvp = 0;

for k = 0:length(waitpdf) - 1
    k_d = min(M, floor(((d / T) - k +delta) / (delta + 2)));
    k_d_CSO = min(M, floor(((d / T) - k+delta-1) / (delta + 2))); % originally it was delta+1
    % Zeta will already be considered whererever this is used. 
    % this is same as k_d = min(M, 1 + floor(((d / T) - 1 - k) / (delta + 1)));
    dvp = dvp + (1-k_d*CSO)*waitpdf(k+1) * prod(pHarq(1:k_d)) + (k_d*CSO)*waitpdf(k + 1) * prod(pHarq(1:k_d_CSO));
end

end
