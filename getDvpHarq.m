function [dvp] = getDvpHarq (d, delta, pHarq, f, T, M, Qmax)
%GETDVPHARQ Function to get the DVP with IR
%   Function description


% k_d = min(M, 1 + floor((d/T - 1) / (delta + 1)));
% dvp_IR_serv = prod(pHarq(1:k_d));

waitpdf = getWaitPmf(pHarq(:), f, M, Qmax);
dvp = 0;
waitpdf
for k = 0:length(waitpdf) - 1
    k_d = min(M, floor(((d / T) - k + delta) / (delta + 1)));
    % Zeta will already be considered whererever this is used. 
    % this is same as k_d = min(M, 1 + floor(((d / T) - 1 - k) / (delta + 1)));
    
    dvp = dvp + waitpdf(k + 1) * prod(pHarq(1:k_d));
end


end
