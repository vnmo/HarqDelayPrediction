function [QLenPmf, P] = getQLenPmfDM(pVector, f, M, Qmax, slotCap)
%GETQLENPMFDM
% Queue length PMF for deterministic-service abstraction with HARQ.
%
% INPUTS
%   pVector  : scalar or length-M vector
%              pVector(m) = failure probability at HARQ attempt m
%   f        : arrival rate (Bernoulli arrival probability per slot)
%   M        : maximum HARQ attempts
%   Qmax     : maximum queue length
%   slotCap  : service capacity per slot
%   pktSize  : packet size
%
% OUTPUTS
%   QLenPmf  : PMF of queue length q = 0..Qmax
%   P        : transition probability matrix
%
% MODEL
%   rho = pktSize / slotCap
%
%   If rho >= 1:
%       one transmission attempt takes either floor(rho) or ceil(rho) slots
%       so an HOL packet occupies the server for multiple slots.
%
%   If rho < 1:
%       one slot can serve either floor(1/rho) or ceil(1/rho) packets
%
%   HARQ:
%       - success at attempt m with probability 1-p(m)
%       - failure with probability p(m)
%       - after attempt M, packet departs regardless (success/drop)
%
% NOTE
%   This is a single-server abstraction for queue-length PMF.
%   It does NOT model parallel HARQ processes separately.
    pktSize = 8*108;
    %% parse pVector
    if isscalar(pVector)
        p = repmat(pVector, 1, M);
    elseif numel(pVector) == M
        p = reshape(pVector, 1, []);
    else
        error('pVector must be scalar or length M.');
    end

    if f < 0 || f > 1
        error('f must be an arrival probability in [0,1].');
    end
    if slotCap <= 0 || pktSize <= 0
        error('slotCap and pktSize must be positive.');
    end

    %% service abstraction
    rho = pktSize / slotCap;

    if rho >= 1
        regime = "multi_slot";
        n1 = floor(rho);
        frac = rho - n1;
        if n1 < 1
            n1 = 1;
        end
        if frac < 1e-14
            serviceTimes = n1;
            serviceProbs = 1;
        else
            serviceTimes = [n1, n1+1];
            serviceProbs = [1-frac, frac];
        end
        Smax = max(serviceTimes);
    else
        regime = "multi_packet";
        mu = 1 / rho; % packets per slot
        k1 = floor(mu);
        frac = mu - k1;
        if k1 < 1
            k1 = 1;
        end
        if frac < 1e-14
            slotCaps = k1;
            slotCapProbs = 1;
        else
            slotCaps = [k1, k1+1];
            slotCapProbs = [1-frac, frac];
        end
        Smax = 1;
    end

    %% state
    % state = (q, m, s)
    % q = 0..Qmax
    % m = 1..M   (HARQ attempt of HOL)
    % s = 1..Smax
    %     multi_slot: remaining slots to finish current attempt
    %     multi_packet: dummy = 1 only
    %
    % when q=0, only (m=1,s=1) is meaningful

    L = (Qmax + 1) * M * Smax;
    idx = @(q,m,s) q*(M*Smax) + (m-1)*Smax + s;

    rows_j = cell(L,1);
    rows_v = cell(L,1);

    %% build transitions
    for q = 0:Qmax
        for m = 1:M
            for s = 1:Smax
                s1 = idx(q,m,s);

                % unreachable empty states
                if q == 0 && ~(m == 1 && s == 1)
                    rows_j{s1} = s1;
                    rows_v{s1} = 1;
                    continue;
                end

                transMap = containers.Map('KeyType','double','ValueType','double');

                if regime == "multi_slot"
                    transMap = stepMultiSlot(q, m, s, p, f, M, Qmax, ...
                        serviceTimes, serviceProbs, Smax, idx, transMap);
                else
                    transMap = stepMultiPacket(q, m, p, f, M, Qmax, ...
                        slotCaps, slotCapProbs, idx, transMap);
                end

                rows_j{s1} = cell2mat(keys(transMap));
                rows_v{s1} = cell2mat(values(transMap));
            end
        end
    end

    %% assemble P
    nnzTotal = 0;
    for i = 1:L
        nnzTotal = nnzTotal + numel(rows_j{i});
    end

    ii = zeros(nnzTotal,1);
    jj = zeros(nnzTotal,1);
    vv = zeros(nnzTotal,1);

    ptr = 0;
    for i = 1:L
        nj = numel(rows_j{i});
        if nj == 0
            continue;
        end
        ii(ptr+1:ptr+nj) = i;
        jj(ptr+1:ptr+nj) = rows_j{i};
        vv(ptr+1:ptr+nj) = rows_v{i};
        ptr = ptr + nj;
    end

    P = sparse(ii, jj, vv, L, L);

    %% sanity
    rowSums = full(sum(P,2));
    maxDev = max(abs(rowSums - 1));
    if maxDev > 1e-10
        warning('Transition matrix row sums deviate from 1. Max deviation = %g', maxDev);
    end

    %% steady-state
    if L <= 2000
        [V,D] = eig(full(P.'));
        vec = [];
        for i = 1:size(D,1)
            if abs(D(i,i)-1) < 1e-10
                vec = real(V(:,i));
                break;
            end
        end
        if isempty(vec)
            error('Could not find eigenvector for eigenvalue 1.');
        end
    else
        try
            opts.tol = 1e-10;
            opts.maxit = 5000;
            vec = eigs(P.', 1, 1, opts);
            vec = real(vec);
        catch
            vec = ones(L,1)/L;
            for it = 1:20000
                vecNew = P.' * vec;
                vecNew = vecNew / sum(vecNew);
                if norm(vecNew - vec, 1) < 1e-12
                    vec = vecNew;
                    break;
                end
                vec = vecNew;
            end
        end
    end

    vec = real(vec);
    vec(vec < 0 & abs(vec) < 1e-12) = 0;
    vec = vec / sum(vec);

    %% marginalize over m,s
    QLenPmf = zeros(1, Qmax+1);
    for q = 0:Qmax
        total = 0;
        for m = 1:M
            for s = 1:Smax
                total = total + vec(idx(q,m,s));
            end
        end
        QLenPmf(q+1) = total;
    end
end

%% ============================================================
function mp = stepMultiSlot(q, m, s, p, f, M, Qmax, serviceTimes, serviceProbs, Smax, idx, mp)
% rho >= 1
% one attempt takes multiple slots

    % Bernoulli arrival at end of slot
    arrivals = [0 1];
    arrProbs = [1-f, f];

    if q == 0
        % empty queue: if arrival occurs, start fresh packet
        for ia = 1:2
            a = arrivals(ia);
            pa = arrProbs(ia);

            if a == 0
                s2 = idx(0,1,1);
                mp = addProb(mp, s2, pa);
            else
                for k = 1:numel(serviceTimes)
                    dur = serviceTimes(k);
                    pd = serviceProbs(k);

                    if dur == 1
                        % attempt completes within this slot
                        if m < M %#ok<*NASGU>
                            sSucc = idx(0,1,1);
                            sFail = idx(1,2,1);
                            mp = addProb(mp, sSucc, pa*pd*(1-p(1)));
                            mp = addProb(mp, sFail, pa*pd*p(1));
                        else
                            s2 = idx(0,1,1);
                            mp = addProb(mp, s2, pa*pd);
                        end
                    else
                        s2 = idx(1,1,dur-1);
                        mp = addProb(mp, s2, pa*pd);
                    end
                end
            end
        end
        return;
    end

    % q >= 1
    if s > 1
        % attempt still ongoing, one slot consumed
        for ia = 1:2
            a = arrivals(ia);
            pa = arrProbs(ia);
            q2 = min(Qmax, q + a);
            s2 = idx(q2, m, s-1);
            mp = addProb(mp, s2, pa);
        end
        return;
    end

    % s == 1: current attempt completes in this slot
    for ia = 1:2
        a = arrivals(ia);
        pa = arrProbs(ia);

        % success
        qSucc = q - 1;
        qSucc = min(Qmax, qSucc + a);
        if qSucc == 0
            sSucc = idx(0,1,1);
        else
            % next HOL starts new attempt with random service time
            for k = 1:numel(serviceTimes)
                dur = serviceTimes(k);
                pd = serviceProbs(k);
                if dur == 1
                    % immediately completes next HOL's first attempt in same slot is NOT allowed
                    sSucc2 = idx(qSucc,1,1);
                else
                    sSucc2 = idx(qSucc,1,dur-1);
                end
                mp = addProb(mp, sSucc2, pa*(1-p(m))*pd);
            end
            sSucc = [];
        end
        if ~isempty(sSucc)
            mp = addProb(mp, sSucc, pa*(1-p(m)));
        end

        % failure
        if m < M
            qFail = min(Qmax, q + a);
            for k = 1:numel(serviceTimes)
                dur = serviceTimes(k);
                pd = serviceProbs(k);
                if dur == 1
                    sFail = idx(qFail, m+1, 1);
                else
                    sFail = idx(qFail, m+1, dur-1);
                end
                mp = addProb(mp, sFail, pa*p(m)*pd);
            end
        else
            % drop after M
            qDone = q - 1;
            qDone = min(Qmax, qDone + a);
            if qDone == 0
                sDone = idx(0,1,1);
                mp = addProb(mp, sDone, pa);
            else
                for k = 1:numel(serviceTimes)
                    dur = serviceTimes(k);
                    pd = serviceProbs(k);
                    if dur == 1
                        sDone2 = idx(qDone,1,1);
                    else
                        sDone2 = idx(qDone,1,dur-1);
                    end
                    mp = addProb(mp, sDone2, pa*pd);
                end
            end
        end
    end
end

%% ============================================================
function mp = stepMultiPacket(q, m, p, f, M, Qmax, slotCaps, slotCapProbs, idx, mp)
% rho < 1
% one slot can serve multiple packets, but same packet retransmission
% within the same slot is NOT allowed.
%
% state only needs s=1 in this regime.

    arrivals = [0 1];
    arrProbs = [1-f, f];

    for ic = 1:numel(slotCaps)
        K = slotCaps(ic);
        pK = slotCapProbs(ic);

        for ia = 1:2
            a = arrivals(ia);
            pa = arrProbs(ia);

            if q == 0
                qIn = a;
                if qIn == 0
                    s2 = idx(0,1,1);
                    mp = addProb(mp, s2, pK*pa);
                    continue;
                end
                mCur = 1;
            else
                qIn = q + a;
                mCur = m;
            end

            % serve up to K distinct packets once each
            qRemain = qIn;
            doneProbMap = containers.Map('KeyType','char','ValueType','double');
            doneProbMap('0_0') = 1.0; % key = "dep_retx"

            xServe = min(K, qIn);

            for t = 1:xServe
                newMap = containers.Map('KeyType','char','ValueType','double');
                keysOld = keys(doneProbMap);

                for kk = 1:numel(keysOld)
                    key = keysOld{kk};
                    probOld = doneProbMap(key);
                    vals = sscanf(key, '%d_%d');
                    dep = vals(1);
                    retx = vals(2);

                    if t == 1 && q > 0
                        thisM = mCur;
                    else
                        thisM = 1;
                    end

                    if thisM < M
                        keySucc = sprintf('%d_%d', dep+1, retx);
                        keyFail = sprintf('%d_%d', dep, retx+1);
                        newMap = addProbStr(newMap, keySucc, probOld*(1-p(thisM)));
                        newMap = addProbStr(newMap, keyFail, probOld*p(thisM));
                    else
                        keyDone = sprintf('%d_%d', dep+1, retx);
                        newMap = addProbStr(newMap, keyDone, probOld);
                    end
                end

                doneProbMap = newMap;
            end

            finalKeys = keys(doneProbMap);
            for kk = 1:numel(finalKeys)
                key = finalKeys{kk};
                probEnd = doneProbMap(key);
                vals = sscanf(key, '%d_%d');
                dep = vals(1);
                retx = vals(2);

                qOut = qIn - dep + retx;
                qOut = min(Qmax, qOut);

                if qOut == 0
                    s2 = idx(0,1,1);
                else
                    if retx > 0
                        s2 = idx(qOut, 2, 1);
                    else
                        s2 = idx(qOut, 1, 1);
                    end
                end

                mp = addProb(mp, s2, pK*pa*probEnd);
            end
        end
    end
end

%% ============================================================
function mp = addProb(mp, key, val)
    if isKey(mp, key)
        mp(key) = mp(key) + val;
    else
        mp(key) = val;
    end
end

function mp = addProbStr(mp, key, val)
    if isKey(mp, key)
        mp(key) = mp(key) + val;
    else
        mp(key) = val;
    end
end