function [waitFunction, pDrop] = getWaitFunctionCSO(pVector, f, M, Qmax, CSO)
% Inputs: conditional PER vector pVector, arrival probability f, attempt limit M, maximum aggregate packet count Qmax, and CSO.
% Outputs: Eq. (11) normalized aggregate waiting-delay function,
% and arrival-discard probability pDrop = pi(Qmax).
% Conditional on admitted q, U is a sum of q independent full-attempt
% workloads. Its CS starting phase is uniform and independent of (q,U).
% This construction is not the DTMC's conditional residual workload, nor
% the physical waiting delay produced by simulate_queue_IR_cso.
% Waiting and service phases are averaged separately in the aggregate model.

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
        'The single-blockage waiting-delay model requires M*Qmax < 1/CSO.');
end

packetCountDistribution = getQLenPmfCSO(p, f, M, Qmax, CSO);
if any(~isfinite(packetCountDistribution)) ...
        || any(packetCountDistribution < -1e-12) ...
        || abs(sum(packetCountDistribution) - 1) > 1e-10
    error('CSO:InvalidQueuePMF', ...
        'The packet-count distribution is invalid.');
end
packetCountDistribution = max(packetCountDistribution, 0);
packetCountDistribution = packetCountDistribution / sum(packetCountDistribution);
pDrop = packetCountDistribution(end);
admissionMass = sum(packetCountDistribution(1:end-1));
if admissionMass <= 0
    error('CSO:NoAdmittedPackets', ...
        'The waiting-delay function is undefined without admissions.');
end

% Algorithm 1: P(L=t)=(1-p(t))*prod(p(1:t-1)) for t<M;
% P(L=M)=prod(p(1:M-1)), including success AND discard at the final attempt.
reachAttempt = [1, cumprod(p(1:M-1))];
attemptPmf = reachAttempt .* [1-p(1:M-1), 1];
underlyingWaitPmf = zeros(1, M * (Qmax-1) + 1);
conditionalWaitPmf = 1;
for q = 0:Qmax-1
    underlyingWaitPmf(1:numel(conditionalWaitPmf)) = ...
        underlyingWaitPmf(1:numel(conditionalWaitPmf)) ...
        + (packetCountDistribution(q+1) / admissionMass) * conditionalWaitPmf;
    if q < Qmax-1
        conditionalWaitPmf = conv(conditionalWaitPmf, [0, attemptPmf]);
    end
end

% Mixing over admitted q can precede the CS correction because both
% operations are linear and beta(k) does not depend on q.
% Eqs. (10)/(11): (1-beta_k)*u_k + beta_(k-1)*u_(k-1).
% Each workload u contributes beta_u to waiting slot u+1. The function
% sums to one by construction; no post-hoc normalization is applied.
k = 0:numel(underlyingWaitPmf);
beta = min(1, max(0, k * CSO));
waitFunction = [(1-beta(1:end-1)) .* underlyingWaitPmf, 0] ...
    + [0, beta(1:end-1) .* underlyingWaitPmf];
end
