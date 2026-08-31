function tests = testCsoDelays
%TESTCSODELAYS Formula and boundary tests; no queueing simulation is run.
tests = functiontests(localfunctions);
end

function testWaitingMatchesEnumeratedCounts(testCase)
p = [0.6, 0.3, 0.1];
M = 3;
Qmax = 4;
CSO = 1/80;
f = 0.4;
qPmf = getQLenPmfCSO(p, f, M, Qmax, CSO);
[actual, drop] = getWaitPmfCSO(p, f, M, Qmax, CSO);
expected = zeros(1, M*(Qmax-1)+2);
for q = 0:Qmax-1
    underlying = enumerateAttempts(q, M, p);
    for k = 0:numel(underlying)-1
        weight = qPmf(q+1) / (1-qPmf(end)) * underlying(k+1);
        expected(k+1) = expected(k+1) + (1-k*CSO)*weight;
        expected(k+2) = expected(k+2) + k*CSO*weight;
    end
end
verifyEqual(testCase, actual, expected, 'AbsTol', 1e-12);
verifyEqual(testCase, drop, qPmf(end), 'AbsTol', 1e-12);
end

function testUnderlyingMatchesOriginalRecursion(testCase)
p = [0.7, 0.3, 0.1, 0.05];
for M = 2:4
    for q = 0:4
        expected = enumerateAttempts(q, M, p);
        actual = arrayfun(@(k) getWaitProbability(k, q, p, M), 0:M*q);
        verifyEqual(testCase, actual, expected, 'AbsTol', 1e-12);
    end
end
end

function testWaitingNormalizationAndSupport(testCase)
for M = [1, 2, 4]
    for Qmax = [1, 4, 16]
        for f = [0, 0.1, 0.6, 0.95]
            for CSO = [0, 1/80]
                [pmf, drop] = getWaitPmfCSO(0.4, f, M, Qmax, CSO);
                verifyEqual(testCase, sum(pmf), 1, 'AbsTol', 1e-12);
                verifyGreaterThanOrEqual(testCase, min(pmf), 0);
                verifyGreaterThanOrEqual(testCase, drop, 0);
                verifyLessThanOrEqual(testCase, drop, 1);
                verifyEqual(testCase, numel(pmf), M*(Qmax-1)+2);
                if f == 0 || Qmax == 1
                    verifyEqual(testCase, pmf(1), 1, 'AbsTol', 1e-12);
                end
            end
        end
    end
end
end

function testZeroIntermediateMassAndFinalShift(testCase)
% Every packet takes all M attempts; positive mass follows initial zeros.
p = [1, 1, 0.3];
M = 3;
Qmax = 4;
CSO = 1/80;
qPmf = getQLenPmfCSO(p, 0.4, M, Qmax, CSO);
pmf = getWaitPmfCSO(p, 0.4, M, Qmax, CSO);
verifyEqual(testCase, pmf(2:3), [0, 0], 'AbsTol', 1e-12);
verifyGreaterThan(testCase, pmf(M+1), 0);
lastK = M*(Qmax-1);
verifyEqual(testCase, pmf(end), ...
    lastK*CSO*qPmf(Qmax)/(1-qPmf(end)), 'AbsTol', 1e-12);
verifyGreaterThan(testCase, pmf(end), 0);
verifyEqual(testCase, sum(pmf), 1, 'AbsTol', 1e-12);
end

function testOneAttemptAndPERShapes(testCase)
for per = [0, 0.4, 1]
    [pmf, drop] = getWaitPmfCSO(per, 0.3, 1, 4, 1/80);
    qPmf = getQLenPmfCSO(per, 0.3, 1, 4, 1/80);
    base = qPmf(1:end-1) / (1-drop);
    beta = (0:3)/80;
    expected = [(1-beta).*base, 0] + [0, beta.*base];
    verifyEqual(testCase, pmf, expected, 'AbsTol', 1e-12);
end
p = [0.5, 0.3, 0.1];
verifyEqual(testCase, getWaitPmfCSO(p, 0.3, 3, 4, 1/80), ...
    getWaitPmfCSO(p.', 0.3, 3, 4, 1/80));
verifyEqual(testCase, getWaitPmfCSO(0.5, 0.3, 3, 4, 1/80), ...
    getWaitPmfCSO([0.5, 0.5, 0.5], 0.3, 3, 4, 1/80));
end

function testServiceThresholdsAndDiscard(testCase)
p = [0.4, 0.2, 0.1];
CSO = 1/80;
x = [-Inf, -1, 0, 1.999, 2, 2.999, 3, 5, 6, 8, 9, Inf];
expected = [1, 1, 1, 1, (1-CSO)*0.4+CSO, ...
    (1-CSO)*0.4+CSO, 0.4, ...
    (1-2*CSO)*0.08+2*CSO*0.4, 0.08, ...
    (1-3*CSO)*0.008+3*CSO*0.08, 0.008, 0.008];
verifyEqual(testCase, getServiceDvpCSO(x, 1, 1, p, 3, CSO), ...
    expected, 'AbsTol', 1e-14);
end

function testServiceAgainstCompletionTimes(testCase)
% Count decoding-completion times instead of using the floor formula.
p = [0.65, 0.2, 0.05];
for tau_d = [0, 0.5, 1.2]
    for tau_f = [0, 0.5, 2.1]
        for CSO = [0, 1/80]
            x = [-Inf, -1, 0:0.25:25, Inf];
            completionTimes = (1:3)*(1+ceil(tau_d)+ceil(tau_f))-ceil(tau_f);
            expected = zeros(size(x));
            for n = 1:numel(x)
                noCS = sum(completionTimes <= x(n));
                withCS = sum(completionTimes+1 <= x(n));
                alpha = min(1, noCS*CSO);
                expected(n) = (1-alpha)*prod(p(1:noCS)) ...
                    + alpha*prod(p(1:withCS));
            end
            actual = getServiceDvpCSO(x, tau_d, tau_f, p, 3, CSO);
            verifyEqual(testCase, actual, expected, 'AbsTol', 1e-13);
            verifyLessThanOrEqual(testCase, max(diff(actual)), 1e-13);
        end
    end
end
end

function testSeparateRoundingAndLegacyCall(testCase)
p = [0.4, 0.2, 0.1];
x = 0:0.25:15;
verifyEqual(testCase, getServiceDvpCSO(x, 0.5, 0.5, p, 3, 1/80), ...
    getServiceDvpCSO(x, 1, 1, p, 3, 1/80));
verifyEqual(testCase, getServiceDvpCSO(x, 0.2, 0.7, p, 3, 1/80), ...
    getServiceDvpCSO(x, 1, 1, p, 3, 1/80));
verifyEqual(testCase, getDvpHarqCSO(x, 1, p, 0.3, 1, 3, 4, 1/80), ...
    getDvpHarqCSO(x, 1, p, 0.3, 1, 3, 4, 1/80, 1));
verifyEqual(testCase, getDvpHarqCSO(x, 0.5, p, 0.3, 1, 3, 4, 1/80, 0.5), ...
    getDvpHarqCSO(x, 1, p, 0.3, 1, 3, 4, 1/80, 1));
end

function testTotalIncludesWaitingTailAndAdmissionOnce(testCase)
p = [0.6, 0.3, 0.1];
T = 1e-3;
CSO = 1/80;
[pmf, drop] = getWaitPmfCSO(p, 0.7, 3, 4, CSO);
d = [-Inf, -1, 0, 1, 2, 2.5, 3, 5, 9, 20, Inf]*T;
expected = zeros(size(d));
for n = 1:numel(d)
    admittedViolation = 0;
    for k = 0:numel(pmf)-1
        if k > floor(d(n)/T)
            conditionalViolation = 1;
        else
            times = (1:3)*3-1;
            noCS = sum(times <= d(n)/T-k);
            withCS = sum(times+1 <= d(n)/T-k);
            alpha = noCS*CSO;
            conditionalViolation = (1-alpha)*prod(p(1:noCS)) ...
                + alpha*prod(p(1:withCS));
        end
        admittedViolation = admittedViolation + pmf(k+1)*conditionalViolation;
    end
    expected(n) = drop + (1-drop)*admittedViolation;
end
actual = getDvpHarqCSO(d, 1, p, 0.7, T, 3, 4, CSO, 1);
verifyEqual(testCase, actual, expected, 'AbsTol', 1e-12);
verifyEqual(testCase, actual(1:4), ones(1,4), 'AbsTol', 1e-12);
verifyGreaterThan(testCase, drop, 0);
verifyEqual(testCase, actual(end), drop+(1-drop)*prod(p), 'AbsTol', 1e-12);
end

function testTotalRangesMonotonicityAndShapes(testCase)
deadlines = [-Inf, -1, 0:0.5:90, Inf];
for M = [1, 2, 4]
    for per = [0, 0.4, 1]
        for CSO = [0, 1/80]
            actual = getDvpHarqCSO(deadlines, 1, per, 0.4, 1, M, 8, CSO, 1);
            verifyGreaterThanOrEqual(testCase, min(actual), -1e-12);
            verifyLessThanOrEqual(testCase, max(actual), 1+1e-12);
            verifyLessThanOrEqual(testCase, max(diff(actual)), 1e-12);
            if per == 1
                verifyEqual(testCase, actual, ones(size(actual)), 'AbsTol', 1e-12);
            end
        end
    end
end
x = [0, 2; 3, Inf];
verifySize(testCase, getDvpHarqCSO(x, 1, 0.4, 0.3, 1, 3, 4, 1/80), size(x));
verifySize(testCase, getServiceDvpCSO(x, 1, 1, 0.4, 3, 1/80), size(x));
verifyEmpty(testCase, getServiceDvpCSO([], 1, 1, 0.4, 3, 1/80));
verifyEmpty(testCase, getDvpHarqCSO([], 1, 0.4, 0.3, 1, 3, 4, 1/80));
end

function testTimeUnitsAtSlotBoundaries(testCase)
p = [0.4, 0.2, 0.1];
slots = 0:100;
reference = getDvpHarqCSO(slots, 1, p, 0.3, 1, 3, 4, 1/80);
for T = [1e-3, 0.25e-3, 0.003]
    actual = getDvpHarqCSO(slots*T, 1, p, 0.3, T, 3, 4, 1/80);
    verifyEqual(testCase, actual, reference, 'AbsTol', 1e-12);
end
end

function testExistingPERTableInputs(testCase)
% Exercise the analytical path with the existing stored link-level inputs.
% This is a compatibility check, not a new link-level or queueing simulation.
for sinrDb = [10, 20]
    for mcs = [0:9, 12, 14, 20]
        p = getLenaHarqBler(1, mcs, sinrDb, 1, 4, 2, 816, 1);
        for CSO = [1/80, 1/320]
            for f = 1./[3, 6, 8, 10]
                values = getDvpHarqCSO([2, 3, 5, 6, 8, 10]*1e-3, ...
                    1, p, f, 1e-3, 4, 16, CSO, 1);
                verifyTrue(testCase, all(isfinite(values)));
                verifyGreaterThanOrEqual(testCase, min(values), -1e-12);
                verifyLessThanOrEqual(testCase, max(values), 1+1e-12);
                verifyLessThanOrEqual(testCase, max(diff(values)), 1e-12);
            end
        end
    end
end
end

function testApproximationDomainChecks(testCase)
verifyError(testCase, @() getWaitPmfCSO(0.4, 0.3, 4, 20, 1/80), ...
    'CSO:WaitingSingleBlockage');
verifyError(testCase, @() getWaitPmfCSO(0.4, 0.3, 4, 21, 1/80), ...
    'CSO:WaitingSingleBlockage');
verifyError(testCase, @() getServiceDvpCSO(10, 1, 2, 0.4, 4, 1/16), ...
    'CSO:ServiceSingleBlockage');
verifyError(testCase, @() getWaitPmfCSO([0.4, 0.2], 0.3, 3, 4, 1/80), ...
    'CSO:InvalidPERLength');
verifyError(testCase, @() getServiceDvpCSO(10, 1, 1, [0.4, 0.2], 3, 1/80), ...
    'CSO:InvalidPERLength');
end

function testInvalidProbabilitiesAndDelays(testCase)
calls = {@() getWaitPmfCSO(NaN, 0.3, 3, 4, 1/80), ...
    @() getWaitPmfCSO(-0.1, 0.3, 3, 4, 1/80), ...
    @() getWaitPmfCSO(0.4, 1.1, 3, 4, 1/80), ...
    @() getServiceDvpCSO(10, -0.5, 1, 0.4, 3, 1/80), ...
    @() getServiceDvpCSO(NaN, 1, 1, 0.4, 3, 1/80), ...
    @() getDvpHarqCSO(10, 1, 0.4, 0.3, 0, 3, 4, 1/80)};
for n = 1:numel(calls)
    threw = false;
    try
        calls{n}();
    catch
        threw = true;
    end
    verifyTrue(testCase, threw);
end
end

function pmf = enumerateAttempts(q, M, p)
% Enumerate all per-packet attempt counts for a small reference case.
pmf = zeros(1, M*q+1);
if q == 0
    pmf(1) = 1;
    return;
end
previous = enumerateAttempts(q-1, M, p);
for count = 1:M
    probability = prod(p(1:count-1));
    if count < M
        probability = probability*(1-p(count));
    end
    for k = 0:numel(previous)-1
        pmf(k+count+1) = pmf(k+count+1)+previous(k+1)*probability;
    end
end
end
