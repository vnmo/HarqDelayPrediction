function tests = testCsoQueueDtmc
tests = functiontests(localfunctions);
end

function testMatrixMatchesTableRows(testCase)
for M = [1,2,4]
    for Qmax = [1,3]
        for f = [0,0.3,1]
            for CSO = [0,1,1/3,1/8]
                p = linspace(0.7,0.2,M);
                [~,actual] = getQLenPmfCSO(p,f,M,Qmax,CSO);
                expected = tableMatrix(p,f,M,Qmax,CSO);
                verifyEqual(testCase,full(actual.P),full(expected),'AbsTol',1e-14);
                verifyEqual(testCase,full(sum(actual.P,2)),ones(size(actual.P,1),1),'AbsTol',1e-14);
                verifyGreaterThanOrEqual(testCase,min(nonzeros(actual.P)),0);
            end
        end
    end
end
end

function testPhaseAdvanceAndUnusedStates(testCase)
M = 3;
Qmax = 4;
xi = 8;
L = (Qmax+1)*M;
[~,model] = getQLenPmfCSO([0.6,0.3,0.1],0.4,M,Qmax,1/xi);
[source,destination] = find(model.P);
sourcePhase = floor((source-1)/L);
destinationPhase = floor((destination-1)/L);
verifyEqual(testCase,destinationPhase,mod(sourcePhase+1,xi));
verifyNumElements(testCase,model.validIndices,(Qmax*M+1)*xi);
for phase = 0:xi-1
    for m = 2:M
        row = phase*L+m;
        next = mod(phase+1,xi)*L+1;
        verifyEqual(testCase,full(model.P(row,next)),1);
        verifyEqual(testCase,nnz(model.P(row,:)),1);
        verifyEqual(testCase,model.stationary(row),0);
    end
end
end

function testFullBufferRejectsReplacingArrival(testCase)
p = [0.6,0.2];
Qmax = 3;
M = 2;
for f = [0,0.3,1]
    [~,model] = getQLenPmfCSO(p,f,M,Qmax,1/4);
    expected = zeros(1,(Qmax+1)*M);
    expected((Qmax-1)*M+1) = 1-p(1);
    expected(Qmax*M+2) = p(1);
    verifyEqual(testCase,full(model.dataSlot(Qmax*M+1,:)),expected,'AbsTol',1e-14);
    expected(:) = 0;
    expected((Qmax-1)*M+1) = 1;
    verifyEqual(testCase,full(model.dataSlot(Qmax*M+M,:)),expected);
end
end

function testOneAttemptManualDistribution(testCase)
for xi = [2,3,8,80]
    for f = [0,0.3,1]
        for p = [0,0.4,1]
            [actual,model] = getQLenPmfCSO(p,f,1,1,1/xi);
            verifyEqual(testCase,actual,[1-f/xi,f/xi],'AbsTol',1e-13);
            verifyEqual(testCase,full(model.dataSlot(1,:)),[1,0]);
            verifyEqual(testCase,full(model.dataSlot(2,:)),[1,0]);
            verifyEqual(testCase,size(model.P),[2*xi,2*xi]);
        end
    end
end
end

function testStationaryAgainstFullLinearSystem(testCase)
for M = [1,2,4]
    for xi = [3,8]
        Qmax = 3;
        [qPmf,model] = getQLenPmfCSO(0.6,0.4,M,Qmax,1/xi);
        indices = model.validIndices;
        P = model.P(indices,indices);
        A = P.'-speye(numel(indices));
        A(end,:) = 1;
        b = zeros(numel(indices),1);
        b(end) = 1;
        expected = A\b;
        verifyEqual(testCase,model.stationary(indices),expected,'AbsTol',2e-12);
        verifyLessThan(testCase,model.stationaryResidual,1e-12);
        verifyEqual(testCase,sum(model.stationary),1,'AbsTol',1e-12);
        verifyEqual(testCase,sum(model.phaseConditional,1),ones(1,xi),'AbsTol',1e-12);
        grouped = reshape(model.stationary,M,Qmax+1,xi);
        marginal = reshape(sum(sum(grouped,1),3),1,[]);
        verifyEqual(testCase,qPmf,marginal,'AbsTol',1e-13);
    end
end
end

function testPhaseChainIsNotAveragedTransitionMatrix(testCase)
M = 2;
Qmax = 3;
xi = 4;
[phasePmf,model] = getQLenPmfCSO([0.6,0.2],0.4,M,Qmax,1/xi);
indices = [1,M+1:(Qmax+1)*M];
averaged = (model.blockedSlot+(xi-1)*model.dataSlot)/xi;
P = averaged(indices,indices);
A = P.'-speye(numel(indices));
A(end,:) = 1;
b = zeros(numel(indices),1);
b(end) = 1;
v = zeros((Qmax+1)*M,1);
v(indices) = A\b;
averagedPmf = sum(reshape(v,M,Qmax+1),1);
verifyGreaterThan(testCase,norm(phasePmf-averagedPmf,Inf),1e-6);
end

function testEndpointAndNoBlockageCases(testCase)
for M = [1,2,4]
    for f = [0,1]
        for p = [0,1]
            for CSO = [0,1,1/2,1/3]
                [actual,model] = getQLenPmfCSO(p,f,M,3,CSO);
                verifyEqual(testCase,sum(actual),1,'AbsTol',1e-12);
                verifyGreaterThanOrEqual(testCase,min(actual),0);
                verifyLessThan(testCase,model.stationaryResidual,1e-12);
                if f == 0
                    verifyEqual(testCase,actual,[1,0,0,0],'AbsTol',1e-12);
                elseif CSO == 1
                    verifyEqual(testCase,actual,[0,0,0,1],'AbsTol',1e-12);
                end
            end
        end
    end
end
verifyEqual(testCase,getQLenPmfCSO(0,1,1,4,0),[1,0,0,0,0]);
end

function testProductionPeriodsAreSparse(testCase)
for xi = [80,320]
    [actual,model] = getQLenPmfCSO([0.4,0.2,0.1,0.05],1/3,4,16,1/xi);
    verifyTrue(testCase,issparse(model.P));
    verifyEqual(testCase,size(model.P),[68*xi,68*xi]);
    verifyLessThanOrEqual(testCase,nnz(model.P),4*68*xi);
    verifyEqual(testCase,sum(actual),1,'AbsTol',1e-12);
    verifyLessThan(testCase,model.stationaryResidual,1e-12);
end
end

function testInputValidationAndShapes(testCase)
p = [0.6,0.3];
verifyEqual(testCase,getQLenPmfCSO(p,0.3,2,3,1/8), ...
    getQLenPmfCSO(p.',0.3,2,3,1/8),'AbsTol',1e-14);
verifyError(testCase,@() getQLenPmfCSO(0.5,0.3,2,3,0.3),'CSO:NonintegerPeriod');
verifyError(testCase,@() getQLenPmfCSO([0.5,0.3],0.3,3,3,1/8),'CSO:InvalidPERLength');
end

function P = tableMatrix(p,f,M,Qmax,CSO)
if CSO == 0
    xi = 1;
else
    xi = round(1/CSO);
end
L = (Qmax+1)*M;
P = sparse(L*xi,L*xi);
for phase = 0:xi-1
    I = double(CSO==0 || phase~=0);
    B = 1-I;
    if M == 1
        add(0,1,0,1,1-B*f);
    else
        add(0,1,0,1,1-B*f-I*f*p(1));
        add(0,1,1,2,I*f*p(1));
    end
    add(0,1,1,1,B*f);
    for m = 2:M
        add(0,m,0,1,1);
    end
    for q = 1:Qmax-1
        if M > 1
            add(q,1,q,1,B*(1-f)+I*f*(1-p(1)));
        end
        for m = 2:M-1
            add(q,m,q,m,B*(1-f));
            add(q,m,q,1,I*f*(1-p(m)));
        end
        for m = 1:M-1
            add(q,m,q+1,m,B*f);
            add(q,m,q,m+1,I*(1-f)*p(m));
            add(q,m,q+1,m+1,I*f*p(m));
            add(q,m,q-1,1,I*(1-f)*(1-p(m)));
        end
        add(q,M,q,1,I*f);
        add(q,M,q-1,1,I*(1-f));
        add(q,M,q,M,B*(1-f));
        add(q,M,q+1,M,B*f);
    end
    for m = 1:M
        add(Qmax,m,Qmax,m,B);
    end
    for m = 1:M-1
        add(Qmax,m,Qmax,m+1,I*p(m));
        add(Qmax,m,Qmax-1,1,I*(1-p(m)));
    end
    add(Qmax,M,Qmax-1,1,I);
end
    function add(q,m,nextQ,nextM,probability)
        row = phase*L+q*M+m;
        col = mod(phase+1,xi)*L+nextQ*M+nextM;
        P(row,col) = P(row,col)+probability;
    end
end
