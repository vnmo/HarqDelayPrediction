function serviceDvp = getServiceDvpCSO(residualSlots, tau_d, tau_f, pHarq, M, CSO)
% Inputs: residual slot budget, decoding delay tau_d, feedback delay tau_f, conditional PER pHarq, attempt limit M, and CSO.
% Output: service-delay upper bound serviceDvp.

validateattributes(residualSlots, {'numeric'}, {'real','nonnan'});
validateattributes(tau_d, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(tau_f, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(M, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(CSO, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(pHarq, {'numeric'}, {'vector','nonempty','real','finite','>=',0,'<=',1});
if isscalar(pHarq)
    pHarq = repmat(pHarq, 1, M);
elseif numel(pHarq) < M
    error('CSO:InvalidPERLength', 'pHarq must contain at least M elements.');
end
p = reshape(pHarq(1:M), 1, []);

tau_d = ceil(tau_d);
tau_f = ceil(tau_f);
RTT = 1 + tau_d + tau_f;
if CSO > 0 && M * RTT >= 1 / CSO
    error('CSO:ServiceSingleBlockage', ...
        'The service-delay upper bound requires M*RTT < 1/CSO.');
end

M_noCS = min(M, max(0, floor((residualSlots + tau_f) / RTT)));
M_CS = min(M, max(0, floor((residualSlots + tau_f - 1) / RTT)));
alpha = min(1, M_noCS * CSO);

failureProduct = [1, cumprod(p)];
failureNoCS = reshape(failureProduct(M_noCS(:)+1), size(residualSlots));
failureCS = reshape(failureProduct(M_CS(:)+1), size(residualSlots));
serviceDvp = (1-alpha) .* failureNoCS + alpha .* failureCS;
end
