function [QLenPmf] = getQLenPmfCSO(pVector, f, M, Qmax, CSO)
%GETQLENPMFCSO Proposed model: auxiliary queue-length PMF with IR HARQ.
% This function constructs the analytical DTMC; it is not a queueing
% simulation. getWaitPmfCSO uses its output in the waiting approximation.
%
% State (q,m): q is the auxiliary maintained-packet count, and m selects
% the conditional PER p(m). It is not an individual HARQ-process index.
% The matrix index is s=q*M+m, with q=0,...,Qmax and m=1,...,M.
% P(s1,s2) gives the probability from source s1 to destination s2.
% CSO=1/xi is the blocked-slot probability in these transitions; CSP phase
% is not a state dimension here. Periodic slot handling is implemented
% separately in simulate_queue_IR_cso using mod(slotId,k_xi).
% Individual decoding/feedback timers are not stored in this DTMC.
%
% pVector : conditional PER vector
% f       : Bernoulli packet-arrival probability
% M       : maximum number of transmission attempts
% Qmax    : maximum auxiliary packet count
% CSO     : control-signaling overhead probability
% QLenPmf(q+1) : stationary probability after summing over m

%% PER vector
if length(pVector) == 1
    p = repmat(pVector, 1, M);
elseif length(pVector) >= M
    p = pVector(1:M);
else
    error("pVector must contain at least M elements.");
end

p = p(:).';

%% Basic parameters
L = (Qmax + 1) * M;

p_prime = 1 - p;
f_prime = 1 - f;
CSO_prime = 1 - CSO;

P = zeros(L, L);

%% Unreachable Cartesian entries: (q=0,m>1)
% These placeholder entries are excluded from the stationary solve below.
for i = 2:M
    P(i, i) = 1;
end

%% Transitions from (q=0,m=1)
q = 0;
m = 1;
s1 = q * M + m;

% (0,1) -> (0,1): no arrival, or an arrival served successfully this slot
P(s1, s1) = ...
    1 - CSO * f - CSO_prime * f * p(1);

% (0,1) -> (1,1): arrival during a CS-blocked slot
q = 1;
m = 1;
s2 = q * M + m;
P(s1, s2) = CSO * f;

% (0,1) -> (1,2): failed first attempt
if M >= 2
    q = 1;
    m = 2;
    s2 = q * M + m;

    P(s1, s2) = CSO_prime * f * p(1);
else
    % When M=1, a failed first attempt is immediately discarded.
    P(s1, s1) = ...
        P(s1, s1) + CSO_prime * f * p(1);
end

%% Transitions for q=1,...,Qmax and m=1,...,M
for q = 1:Qmax
    for m = 1:M

        s1 = q * M + m;

        %% Interior states: q<Qmax and m<M
        if q < Qmax && m < M

            % (q,m) -> (q,m+1): failed attempt, no arrival
            s2 = q * M + (m + 1);
            P(s1, s2) = ...
                CSO_prime * f_prime * p(m);

            % (q,m) -> (q+1,m+1): failed attempt, one arrival
            s2 = (q + 1) * M + (m + 1);
            P(s1, s2) = ...
                CSO_prime * f * p(m);
        end

        %% CS-blocked slot for q<Qmax
        if q < Qmax

            % (q,m) -> (q+1,m): blocked slot, one arrival
            s2 = (q + 1) * M + m;
            P(s1, s2) = CSO * f;

            % (q,m) -> (q,m): blocked slot, no arrival
            s2 = q * M + m;
            P(s1, s2) = ...
                P(s1, s2) + CSO * f_prime;
        end

        %% Successful decoding before the final attempt
        if m < M

            % (q,m) -> (q,1): success with a replacing arrival
            s2 = q * M + 1;
            P(s1, s2) = ...
                P(s1, s2) + ...
                CSO_prime * f * p_prime(m);

            % (q,m) -> (q-1,1): success without an arrival
            s2 = (q - 1) * M + 1;
            P(s1, s2) = ...
                CSO_prime * f_prime * p_prime(m);
        end

        %% Final attempt: success or discard
        if m == M

            % (q,M) -> (q,1): success/discard with a replacing arrival
            s2 = q * M + 1;
            P(s1, s2) = ...
                P(s1, s2) + CSO_prime * f;

            % (q,M) -> (q-1,1): success/discard without an arrival
            s2 = (q - 1) * M + 1;
            P(s1, s2) = ...
                CSO_prime * f_prime;
        end

        %% Full-buffer boundary
        if q == Qmax && m < M

            % (Qmax,m) -> (Qmax,m+1): failure; q cannot increase
            s2 = q * M + (m + 1);
            P(s1, s2) = ...
                CSO_prime * p(m);
        end

        if q == Qmax
            % Full and blocked: keep (Qmax,m), including any arrival.
            P(s1, s1) = ...
                P(s1, s1) + 1 - sum(P(s1, :));
        end
    end
end

%% Solve the stationary distribution over the reachable state space
%
% The full Cartesian indexing is retained, but entries (q=0,m>1)
% are excluded from the stationary-state calculation.

validIndices = [1, (M + 1):L];
Pvalid = P(validIndices, validIndices);

[V, D] = eig(Pvalid.');

unitIndex = find( ...
    abs(diag(D) - 1) < sqrt(eps), ...
    1, ...
    "first");

if isempty(unitIndex)
    error("No unit eigenvalue was found.");
end

vecValid = real(V(:, unitIndex));
vecValid = vecValid / sum(vecValid);

%% Restore the full Cartesian indexing
vec = zeros(L, 1);
vec(validIndices) = vecValid;

%% Marginalize over m
% The full-buffer entry is used as an auxiliary overflow approximation
% by getWaitPmfCSO, not as a measured physical buffer probability.
QLenPmf = sum(reshape(vec, M, Qmax + 1), 1);
QLenPmf = QLenPmf(:).';

end
