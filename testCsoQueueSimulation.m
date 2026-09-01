function tests = testCsoQueueSimulation
tests = functiontests(localfunctions);
end

function setup(testCase)
state = rng;
testCase.addTeardown(@() rng(state));
rng(30,'twister');
end

function testSuccessRetainsCapacityUntilDecoding(testCase)
[q,count,latency,waiting,drops,s] = simulate_queue_IR_cso(true,1,0,1,Inf,3,1,1,true);
verifyEqual(testCase,q,[1,1,1]);
verifyEqual(testCase,count,[1,1]);
verifyEqual(testCase,latency,[2,2]);
verifyEqual(testCase,waiting,[0,0]);
verifyEqual(testCase,drops,1);
verifyEqual(testCase,s.arrivalSlot,[1,3]);
verifyEqual(testCase,s.completionBoundary,[3,5]);
verifyEqual(testCase,s.populationBeforeArrival,[0,1,0]);
verifyEqual(testCase,s.generatedPackets,3);
verifyEqual(testCase,s.unfinishedPackets,0);
verifyEqual(testCase,s.finalPopulation,0);
end

function testFinalFailureRetainsCapacityUntilDecoding(testCase)
[~,count,latency,~,drops,s] = simulate_queue_IR_cso(true,1,1,1,Inf,3,1,1,true);
verifyEqual(testCase,count,[1,1]);
verifyEqual(testCase,latency,[Inf,Inf]);
verifyEqual(testCase,drops,1);
verifyEqual(testCase,s.completionBoundary,[3,5]);
verifyEqual(testCase,s.discardedPackets,2);
verifyEqual(testCase,s.successfulPackets,0);
end

function testPendingRetriesRemainFifoAcrossBlockage(testCase)
[~,count,latency,~,drops,s] = simulate_queue_IR_cso(true,1,[1,0],1,4,2,4,2,true);
verifyEqual(testCase,s.txSlot,[1,2,5,6]);
verifyEqual(testCase,s.txPacketId,[1,2,1,2]);
verifyEqual(testCase,s.txReadySince,[1,2,4,5]);
verifyEqual(testCase,count,[2,2]);
verifyEqual(testCase,latency,[6,6]);
verifyEqual(testCase,drops,0);
verifyEqual(testCase,s.completionBoundary,[7,8]);
end

function testRetriesHavePriorityOverFirstTransmissions(testCase)
[~,~,~,~,~,s] = simulate_queue_IR_cso(true,1,[1,0],1,Inf,5,8,2,true);
atFour = find(s.txSlot==4,1);
verifyEqual(testCase,s.txPacketId(atFour),1);
verifyEqual(testCase,s.txAttempt(atFour),2);
verifyGreaterThan(testCase,s.firstTxSlot(4),4);
end

function testFeedbackAvailableAtBoundaryCanTransmit(testCase)
for decoding = [0,0.5,1.2]
    for feedback = [0,0.5,2.1]
        [~,~,latency,~,~,s] = simulate_queue_IR_cso(true,1,[1,0],feedback,Inf,1,2,2,true,decoding);
        RTT = 1+ceil(decoding)+ceil(feedback);
        verifyEqual(testCase,s.txSlot,[1,1+RTT]);
        verifyEqual(testCase,latency,2*RTT-ceil(feedback));
    end
end
end

function testDecodingContinuesInBlockedSlots(testCase)
[~,~,latency,~,drops,s] = simulate_queue_IR_cso(true,1,0,1,3,3,1,1,true);
verifyEqual(testCase,s.arrivalSlot,[1,3]);
verifyEqual(testCase,s.firstTxSlot,[1,4]);
verifyEqual(testCase,s.completionBoundary,[3,6]);
verifyEqual(testCase,latency,[2,3]);
verifyEqual(testCase,drops,1);
verifyFalse(testCase,any(mod(s.txSlot,3)==0));
end

function testArrivalInFinalBlockedSlotIsDrained(testCase)
[q,count,latency,waiting,drops,s] = simulate_queue_IR_cso(true,0.5,0,1,2,2,2,1,true);
verifyEqual(testCase,q,[0,1]);
verifyEqual(testCase,count,1);
verifyEqual(testCase,latency,3);
verifyEqual(testCase,waiting,1);
verifyEqual(testCase,s.arrivalSlot,2);
verifyEqual(testCase,s.txSlot,3);
verifyEqual(testCase,s.completionBoundary,5);
verifyEqual(testCase,s.drainSlots,3);
verifyEqual(testCase,drops,0);
end

function testAllGeneratedPacketsEnterDvp(testCase)
for p = {[0.8,0.4,0.1],[1,1,1]}
    [q,count,latency,~,drops,s] = simulate_queue_IR_cso(true,0.7,p{1},1,4,300,3,3,false);
    verifyEqual(testCase,numel(latency)+drops,s.generatedPackets);
    verifyEqual(testCase,numel(latency),s.admittedPackets);
    verifyEqual(testCase,sum(isinf(latency)),s.discardedPackets);
    verifyEqual(testCase,s.successfulPackets+s.discardedPackets,s.admittedPackets);
    verifyGreaterThanOrEqual(testCase,min(q),0);
    verifyLessThanOrEqual(testCase,max(q),3);
    verifyGreaterThanOrEqual(testCase,min(count),1);
    verifyLessThanOrEqual(testCase,max(count),3);
    verifyEqual(testCase,s.unfinishedPackets,0);
    verifyEqual(testCase,s.finalPopulation,0);
    verifyLessThanOrEqual(testCase,s.drainSlots,3*3*(3+1)+1);
    if all(p{1}==1)
        verifyEqual(testCase,count,3*ones(size(count)));
        verifyEqual(testCase,(sum(latency>8)+drops)/s.generatedPackets,1);
    end
end
end

function testDelayEndpointExcludesFinalAck(testCase)
[~,~,latency,~,drops,s] = simulate_queue_IR_cso(true,1,0,5,Inf,1,2,1,true);
verifyEqual(testCase,latency,2);
verifyEqual(testCase,(sum(latency>2)+drops)/s.generatedPackets,0);
verifyEqual(testCase,(sum(latency>1)+drops)/s.generatedPackets,1);
end

function testLegacyOutputsAndOptionalDetailsMatch(testCase)
rng(17,'twister');
[q1,c1,l1,w1,d1] = simulate_queue_IR_cso(true,0.3,[0.7,0.2],1,8,100,4,2,false);
rng(17,'twister');
[q2,c2,l2,w2,d2,~] = simulate_queue_IR_cso(true,0.3,[0.7,0.2],1,8,100,4,2,false,1);
verifyEqual(testCase,{q1,c1,l1,w1,d1},{q2,c2,l2,w2,d2});
end

function testEmptyInputsAndInvalidParameters(testCase)
for N = [0,4]
    [q,c,l,w,d,s] = simulate_queue_IR_cso(true,0,0.4,1,8,N,2,1,false);
    verifyEqual(testCase,q,zeros(1,N));
    verifyEmpty(testCase,c);
    verifyEmpty(testCase,l);
    verifyEmpty(testCase,w);
    verifyEqual(testCase,d,0);
    verifyEqual(testCase,s.generatedPackets,0);
end
verifyError(testCase,@() simulate_queue_IR_cso(true,0.4,0.2,1,8,10,2,1,true), ...
    'CSO:DeterministicArrivalPeriod');
verifyError(testCase,@() simulate_queue_IR_cso(true,0.3,0.2,1,1,10,2,1,false), ...
    'CSO:SimulationPeriod');
verifyError(testCase,@() simulate_queue_IR_cso(true,0.3,[0.2,0.1],1,8,10,2,1,false), ...
    'CSO:InvalidPERLength');
end

function testAgainstIndependentPacketLedger(testCase)
for seed = [1,30]
    for priority = [false,true]
        for stages = [0,1,2]
            p = [0.9,0.5,0.2];
            rng(seed,'twister');
            [q,c,l,w,d,s] = simulate_queue_IR_cso(priority,0.6,p,stages,4,150,3,3,false,stages);
            rng(seed,'twister');
            ref = packetLedger(priority,0.6,p,stages,4,150,3,stages);
            verifyEqual(testCase,q,ref.queue);
            verifyEqual(testCase,c,ref.count);
            verifyEqual(testCase,l,ref.latency);
            verifyEqual(testCase,w,ref.waiting);
            verifyEqual(testCase,d,ref.drops);
            verifyEqual(testCase,s.txSlot,ref.txSlot);
            verifyEqual(testCase,s.txPacketId,ref.txId);
            verifyEqual(testCase,s.completionBoundary,ref.finish);
        end
    end
end
end

function result = packetLedger(priority,f,p,feedback,xi,N,Qmax,decoding)
status = zeros(1,N); % Zero denotes an unused packet record.
arrival = zeros(1,N);
count = zeros(1,N);
first = zeros(1,N);
finish = zeros(1,N);
readyAt = zeros(1,N);
decodeAt = zeros(1,N);
outcome = false(1,N);
discarded = false(1,N);
numPackets = 0;
queue = zeros(1,N);
drops = 0;
txSlot = [];
txId = [];
t = 1;
while t <= N || any(status>0 & status<4)
    for id = find(status==2 & decodeAt==t)
        if outcome(id) || count(id)==numel(p)
            status(id) = 4;
            finish(id) = t;
            discarded(id) = ~outcome(id);
        else
            status(id) = 3;
            readyAt(id) = t+feedback;
        end
    end
    status(status==3 & readyAt==t) = 1;
    if t <= N && rand < f
        if sum(status>0 & status<4) == Qmax
            drops = drops+1;
        else
            numPackets = numPackets+1;
            id = numPackets;
            status(id) = 1;
            arrival(id) = t;
            count(id) = 0;
            first(id) = 0;
            finish(id) = 0;
            readyAt(id) = t;
            decodeAt(id) = 0;
            outcome(id) = false;
            discarded(id) = false;
        end
    end
    ready = find(status==1);
    if mod(t,xi)~=0 && ~isempty(ready)
        if priority
            keys = [-double(count(ready)>0).',readyAt(ready).',ready.'];
        else
            keys = [readyAt(ready).',ready.'];
        end
        [~,order] = sortrows(keys);
        id = ready(order(1));
        if count(id)==0
            first(id) = t;
        end
        count(id) = count(id)+1;
        status(id) = 2;
        decodeAt(id) = t+1+decoding;
        outcome(id) = rand < 1-p(count(id));
        txSlot(end+1) = t; %#ok<AGROW>
        txId(end+1) = id; %#ok<AGROW>
    end
    if t <= N
        queue(t) = sum(status>0 & status<4);
    end
    t = t+1;
end
packets = 1:numPackets;
latency = finish(packets)-arrival(packets);
latency(discarded(packets)) = Inf;
result = struct('queue',queue,'count',count(packets),'latency',latency, ...
    'waiting',first(packets)-arrival(packets),'drops',drops, ...
    'txSlot',txSlot,'txId',txId,'finish',finish(packets));
end
