function [probDistr, pDrop] = getWaitPmfCSO(pVector, f, M, Qmax, CSO)
%GETWAITPMFCSO Proposed model: approximate waiting PMF given admission.
% This analytical function uses getQLenPmfCSO's auxiliary (q,m) DTMC;
% it does not simulate packet queues or advance a CSP-phase clock.
% probDistr(k+1) is the probability of waiting k slots, given admission.
% pDrop is pi_Qmax from the existing auxiliary queue model. It is used as
% an overflow approximation, not as a physical buffer-occupancy identity.
% pVector contains conditional PERs; a scalar is repeated for M attempts.
% CSO = 1/xi; CSO = 0 disables blockage. The model requires M*Qmax < xi.

validateattributes(M, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(Qmax, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(f, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(CSO, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(pVector, {'numeric'}, {'vector','nonempty','real','finite','>=',0,'<=',1});
if isscalar(pVector)
    pVector = repmat(pVector, 1, M);
elseif numel(pVector) < M
    error('CSO:InvalidPERLength', 'pVector must contain at least M elements.');
end
p = reshape(pVector(1:M), 1, []);
if CSO > 0 && M * Qmax >= 1 / CSO
    error('CSO:WaitingSingleBlockage', ...
        'The waiting-delay approximation requires M*Qmax < 1/CSO.');
end

% Obtain the Proposed model's auxiliary queue PMF; remove roundoff negatives.
qLenPmf = getQLenPmfCSO(p, f, M, Qmax, CSO);
if any(~isfinite(qLenPmf)) || any(qLenPmf < -1e-12) ...
        || abs(sum(qLenPmf) - 1) > 1e-10
    error('CSO:InvalidQueuePMF', 'The auxiliary queue PMF is not a probability distribution.');
end
qLenPmf = max(qLenPmf, 0);
qLenPmf = qLenPmf / sum(qLenPmf);
pDrop = qLenPmf(end);
admissionMass = sum(qLenPmf(1:end-1));
if admissionMass <= 0
    error('CSO:NoAdmittedPackets', 'The conditional waiting PMF is undefined without admissions.');
end

% A packet consumes 1,...,M attempts. At M it leaves on success or discard.
% Convolution evaluates the same attempt-count distribution as the waiting
% recursion, including M=1 and zero-probability intermediate counts.
reachAttempt = [1, cumprod(p(1:M-1))];
attemptPmf = reachAttempt .* [1-p(1:M-1), 1];
underlyingProbDistr = zeros(1, M * (Qmax-1) + 1);
conditionalPmf = 1; % q=0: no waiting
for q = 0:Qmax-1
    underlyingProbDistr(1:numel(conditionalPmf)) = ...
        underlyingProbDistr(1:numel(conditionalPmf)) ...
        + (qLenPmf(q+1) / admissionMass) * conditionalPmf;
    if q < Qmax-1
        conditionalPmf = conv(conditionalPmf, [0, attemptPmf]);
    end
end

% fbar(k) = (1-beta_k)*fu(k) + beta_(k-1)*fu(k-1).
% beta is a phase-averaged blockage weight, not a stored CSP-phase state.
% Allocate the extra slot so the final shifted probability is retained.
k = 0:numel(underlyingProbDistr)-1;
beta = min(1, max(0, k * CSO));
probDistr = [(1-beta) .* underlyingProbDistr, 0] ...
    + [0, beta .* underlyingProbDistr];
end
