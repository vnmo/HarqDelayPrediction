function [QLenPmf, details] = getQLenPmfCSO(pVector, f, M, Qmax, CSO)
% Inputs: conditional PER vector pVector, arrival probability f, attempt limit M, maximum aggregate packet count Qmax, and CSO.
% Outputs: packet-count distribution QLenPmf and optional state details.

validateattributes(M, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(Qmax, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(f, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(CSO, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(pVector, {'numeric'}, {'vector','nonempty','real','finite','>=',0,'<=',1});
if isscalar(pVector)
    p = repmat(pVector, 1, M);
elseif numel(pVector) >= M
    p = reshape(pVector(1:M), 1, []);
else
    error('CSO:InvalidPERLength', 'pVector must contain at least M elements.');
end
if CSO == 0
    xi = 1;
else
    rawPeriod = 1 / CSO;
    xi = round(rawPeriod);
    if ~isfinite(rawPeriod) || abs(rawPeriod-xi) > 32*eps(rawPeriod)
        error('CSO:NonintegerPeriod', 'CSO must be zero or the reciprocal of an integer CSP length.');
    end
end

L = (Qmax+1)*M;
dataSlot = zeros(L);
blockedSlot = zeros(L);
for m = 2:M
    dataSlot(m,1) = 1;
    blockedSlot(m,1) = 1;
end
for q = 0:Qmax
    if q == 0
        attempts = 1;
    else
        attempts = 1:M;
    end
    for m = attempts
        row = q*M+m;

        if q < Qmax
            blockedSlot(row,row) = 1-f;
            blockedSlot(row,(q+1)*M+m) = f;
        else
            blockedSlot(row,row) = 1;
        end

        if q == 0
            if M == 1
                dataSlot(row,1) = 1;
            else
                dataSlot(row,1) = 1-f*p(1);
                dataSlot(row,M+2) = f*p(1);
            end
        elseif q == Qmax
            if m == M
                dataSlot(row,(q-1)*M+1) = 1;
            else
                dataSlot(row,(q-1)*M+1) = 1-p(m);
                dataSlot(row,q*M+m+1) = p(m);
            end
        else
            for arrival = 0:1
                arrivalProbability = (1-arrival)*(1-f)+arrival*f;
                if m == M
                    next = (q-1+arrival)*M+1;
                    dataSlot(row,next) = dataSlot(row,next)+arrivalProbability;
                else
                    next = (q-1+arrival)*M+1;
                    dataSlot(row,next) = dataSlot(row,next)+arrivalProbability*(1-p(m));
                    next = (q+arrival)*M+m+1;
                    dataSlot(row,next) = dataSlot(row,next)+arrivalProbability*p(m);
                end
            end
        end
    end
end

validBase = [1, M+1:L];
D = dataSlot(validBase,validBase);
if CSO == 0
    firstSlot = D;
else
    firstSlot = blockedSlot(validBase,validBase);
end

cycle = firstSlot*(D^(xi-1));
phaseZero = stationaryFromEmpty(cycle);
phaseConditional = zeros(numel(validBase),xi);
phaseConditional(:,1) = phaseZero;
nextPhase = firstSlot.'*phaseZero;
for phase = 2:xi
    phaseConditional(:,phase) = nextPhase;
    nextPhase = D.'*nextPhase;
end
cycleResidual = norm(nextPhase-phaseZero,Inf);
if cycleResidual > 1e-9
    error('CSO:StationaryResidual', 'The CSP-cycle stationary residual exceeds tolerance.');
end
basePmf = zeros(L,1);
basePmf(validBase) = sum(phaseConditional,2)/xi;
QLenPmf = sum(reshape(basePmf,M,Qmax+1),1);

if nargout > 1
    phaseAdvance = sparse(1:xi,[2:xi,1],ones(1,xi),xi,xi);
    P = kron(phaseAdvance,sparse(dataSlot));
    if CSO > 0
        P(1:L,:) = 0;
        nextColumns = mod(1,xi)*L+(1:L);
        P(1:L,nextColumns) = sparse(blockedSlot);
    end
    validIndices = reshape(validBase(:)+(0:xi-1)*L,[],1);
    stationary = zeros(L*xi,1);
    stationary(validIndices) = phaseConditional(:)/xi;
    details = struct('P',P,'stationary',stationary, ...
        'validIndices',validIndices,'xi',xi,'phaseConditional',phaseConditional, ...
        'dataSlot',sparse(dataSlot),'blockedSlot',sparse(blockedSlot), ...
        'stationaryResidual',norm(P.'*stationary-stationary,Inf));
end
end

function stationary = stationaryFromEmpty(P)
P = P./sum(P,2);
graphP = digraph(sparse(P>0));
reachable = sort(bfsearch(graphP,1));
R = P(reachable,reachable);
classId = conncomp(digraph(sparse(R>0)),'Type','strong');
closed = {};
for id = unique(classId)
    members = find(classId==id);
    if ~any(R(members,classId~=id)>0,'all')
        closed{end+1} = members; %#ok<AGROW>
    end
end
if isempty(closed)
    error('CSO:NoClosedClass', 'No closed communicating class was found.');
end
transient = true(numel(reachable),1);
classPmf = cell(size(closed));
for j = 1:numel(closed)
    members = closed{j};
    transient(members) = false;
    A = R(members,members).'-eye(numel(members));
    A(end,:) = 1;
    b = zeros(numel(members),1);
    b(end) = 1;
    classPmf{j} = normalizeProbability(A\b);
end
start = find(reachable==1,1);
weights = zeros(1,numel(closed));
if transient(start)
    transientIndices = find(transient);
    exits = zeros(sum(transient),numel(closed));
    for j = 1:numel(closed)
        exits(:,j) = sum(R(transient,closed{j}),2);
    end
    absorption = (eye(sum(transient))-R(transient,transient))\exits;
    weights = normalizeProbability(absorption(transientIndices==start,:).').';
else
    for j = 1:numel(closed)
        weights(j) = any(closed{j}==start);
    end
end
stationary = zeros(size(P,1),1);
for j = 1:numel(closed)
    stationary(reachable(closed{j})) = weights(j)*classPmf{j};
end
end

function p = normalizeProbability(p)
if any(~isfinite(p)) || any(p < -1e-10) || sum(p) <= 0
    error('CSO:InvalidStationaryPMF', 'The stationary solution is not a probability vector.');
end
p = max(p,0);
p = p/sum(p);
end
