function tests = testManuscriptAggregate
tests = functiontests(localfunctions);
end

function testWaitingFunctionAgainstEnumeratedCSPhases(testCase)
% Enumerate periodic blocked slots directly, independently of beta mixing.
for M = 1:3
    for Qmax = [1,4]
        for f = [0,0.4,0.9]
            xi = 20;
            p = [0.7,0.3,0.1];
            piq = getQLenPmfCSO(p,f,M,Qmax,1/xi);
            [actual,drop] = getWaitFunctionCSO(p,f,M,Qmax,1/xi);
            exact = zeros(size(actual));
            for q = 0:Qmax-1
                weight = piq(q+1)/(1-drop);
                for u = 0:M*q
                    mass = getWaitProbability(u,q,p,M)*weight;
                    for phase = 0:xi-1
                        sent = 0; elapsed = 0;
                        while sent < u
                            if mod(phase+elapsed,xi) ~= 0
                                sent = sent+1;
                            end
                            elapsed = elapsed+1;
                        end
                        exact(elapsed+1) = exact(elapsed+1)+mass/xi;
                    end
                end
            end
            verifyEqual(testCase,sum(exact),1,'AbsTol',2e-12);
            verifyEqual(testCase,actual,exact,'AbsTol',2e-12);
            verifyEqual(testCase,sum(actual),1,'AbsTol',2e-12);
        end
    end
end
end

function testTheoremAgainstCompletionTimeCounting(testCase)
% Independently sum Eq. (11) and count feasible decoding completions.
p = [0.6,0.3,0.1,0.05];
for M = [1,2,4]
    for Qmax = [1,4]
        for tau = [0.5,1,1.2]
            xi = 80; f = 0.7; T = 1e-3;
            slots = [-1,0:0.5:20,Inf];
            piq = getQLenPmfCSO(p,f,M,Qmax,1/xi);
            wait = zeros(1,M*(Qmax-1)+2);
            for k = 0:numel(wait)-1
                for q = 0:Qmax-1
                    uk = getWaitProbability(k,q,p,M);
                    previous = getWaitProbability(k-1,q,p,M);
                    wait(k+1) = wait(k+1)+piq(q+1)/(1-piq(end)) ...
                        *((1-k/xi)*uk+(max(0,k-1)/xi)*previous);
                end
            end
            expected = zeros(size(slots));
            done = (1:M)*(1+ceil(tau)+ceil(tau))-ceil(tau);
            for n = 1:numel(slots)
                admitted = 0;
                for k = 0:numel(wait)-1
                    if k>floor(slots(n))
                        value = 1;
                    else
                        noCS = sum(done<=slots(n)-k);
                        withCS = sum(done+1<=slots(n)-k);
                        value = (1-noCS/xi)*prod(p(1:noCS)) ...
                            +(noCS/xi)*prod(p(1:withCS));
                    end
                    admitted = admitted+wait(k+1)*value;
                end
                expected(n) = min(1,piq(end)+(1-piq(end))*admitted);
            end
            actual = getDvpHarqCSO(slots*T,tau,p,f,T,M,Qmax,1/xi,tau);
            verifyEqual(testCase,actual,expected,'AbsTol',2e-12);
        end
    end
end
end
