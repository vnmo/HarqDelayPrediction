function [queueSize, TxCount, totalLatencyInSlotsVector, waitLatencyInslotsVector, NumDroppedPkts, details] = ...
    simulate_queue_IR_cso( ...
    isReTxPrioritised, p_arrival, perVector, feedbackDelay, k_xi, ...
    MaxNumSlots, maxQueueSize, MaxTxPerPkt, isDetArr, decodingDelay)
% Inputs: retry priority, arrival probability, conditional PER vector, feedback delay, CSP length, slot count, maximum aggregate packet count, attempt limit, arrival mode, and decoding delay.
% Outputs: aggregate packet count, attempt counts, total and accumulated waiting delays, arrival-discard count, and event details.

if nargin < 10
    decodingDelay = 1;
end
validateattributes(isReTxPrioritised, {'logical','numeric'}, {'scalar','binary'});
validateattributes(isDetArr, {'logical','numeric'}, {'scalar','binary'});
validateattributes(p_arrival, {'numeric'}, {'scalar','real','finite','>=',0,'<=',1});
validateattributes(MaxNumSlots, {'numeric'}, {'scalar','real','finite','integer','nonnegative'});
validateattributes(maxQueueSize, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(MaxTxPerPkt, {'numeric'}, {'scalar','real','finite','integer','positive'});
validateattributes(feedbackDelay, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(decodingDelay, {'numeric'}, {'scalar','real','finite','nonnegative'});
validateattributes(k_xi, {'numeric'}, {'scalar','real','positive'});
if ~isinf(k_xi) && (k_xi ~= round(k_xi) || k_xi < 2)
    error('CSO:SimulationPeriod', 'The CSP length must be an integer >=2, or Inf without blockage.');
end
validateattributes(perVector, {'numeric'}, {'vector','nonempty','real','finite','>=',0,'<=',1});
if numel(perVector) ~= MaxTxPerPkt
    error('CSO:InvalidPERLength', 'perVector must contain one conditional PER per allowed attempt.');
end
p_ir = reshape(perVector,1,[]);
feedbackDelay = ceil(feedbackDelay);
decodingDelay = ceil(decodingDelay);
RTT = 1+decodingDelay+feedbackDelay;
if isDetArr && p_arrival > 0
    arrivalPeriod = 1/p_arrival;
    if abs(arrivalPeriod-round(arrivalPeriod)) > 32*eps(arrivalPeriod)
        error('CSO:DeterministicArrivalPeriod', 'Deterministic arrivals require an integer cycle length.');
    end
    arrivalPeriod = round(arrivalPeriod);
else
    arrivalPeriod = Inf;
end

queueSize = zeros(1,MaxNumSlots);
arrSlot = zeros(1,MaxNumSlots);
TxCount = zeros(1,MaxNumSlots);
firstTxSlot = zeros(1,MaxNumSlots);
completionBoundary = zeros(1,MaxNumSlots);
discarded = false(1,MaxNumSlots);
accumulatedWait = zeros(1,MaxNumSlots);
readySince = zeros(1,MaxNumSlots);
NumDroppedPkts = 0;
numGenerated = 0;
numAdmitted = 0;
numPacketsInSystem = 0;
readyIds = [];
numReadyRetries = 0;

decodeRingLength = decodingDelay+2;
decodePacket = zeros(1,decodeRingLength);
decodeSuccess = false(1,decodeRingLength);
feedbackRingLength = feedbackDelay+1;
feedbackPacket = zeros(1,feedbackRingLength);

recordDetails = nargout > 5;
if recordDetails
    populationBeforeArrival = zeros(1,MaxNumSlots);
    txSlotLog = [];
    txPacketLog = [];
    txAttemptLog = [];
    txReadySinceLog = [];
end

slotId = 1;
drainLimit = MaxNumSlots+maxQueueSize*MaxTxPerPkt*(RTT+1)+1;
while slotId <= MaxNumSlots || numPacketsInSystem > 0
    if slotId > drainLimit
        error('CSO:DrainDidNotFinish', 'The admitted population did not drain within the event bound.');
    end

    decodeIndex = mod(slotId-1,decodeRingLength)+1;
    pktId = decodePacket(decodeIndex);
    if pktId ~= 0
        succeeded = decodeSuccess(decodeIndex);
        decodePacket(decodeIndex) = 0;
        if succeeded || TxCount(pktId) == MaxTxPerPkt
            completionBoundary(pktId) = slotId;
            discarded(pktId) = ~succeeded;
            numPacketsInSystem = numPacketsInSystem-1;
        else
            readyBoundary = slotId+feedbackDelay;
            feedbackIndex = mod(readyBoundary-1,feedbackRingLength)+1;
            if feedbackPacket(feedbackIndex) ~= 0
                error('CSO:FeedbackEventCollision', 'Two packets were scheduled for the same feedback entry.');
            end
            feedbackPacket(feedbackIndex) = pktId;
        end
    end

    feedbackIndex = mod(slotId-1,feedbackRingLength)+1;
    pktId = feedbackPacket(feedbackIndex);
    if pktId ~= 0
        feedbackPacket(feedbackIndex) = 0;
        if isReTxPrioritised
            readyIds = [readyIds(1:numReadyRetries), pktId, readyIds(numReadyRetries+1:end)];
            numReadyRetries = numReadyRetries+1;
        else
            readyIds(end+1) = pktId; %#ok<AGROW>
        end
        readySince(pktId) = slotId;
    end

    if slotId <= MaxNumSlots
        if recordDetails
            populationBeforeArrival(slotId) = numPacketsInSystem;
        end
        if p_arrival == 0
            isArrival = false;
        elseif isDetArr
            isArrival = mod(slotId,arrivalPeriod)==0;
        else
            isArrival = rand < p_arrival;
        end
        if isArrival
            numGenerated = numGenerated+1;
            if numPacketsInSystem >= maxQueueSize
                NumDroppedPkts = NumDroppedPkts+1;
            else
                numAdmitted = numAdmitted+1;
                numPacketsInSystem = numPacketsInSystem+1;
                arrSlot(numAdmitted) = slotId;
                readyIds(end+1) = numAdmitted; %#ok<AGROW>
                readySince(numAdmitted) = slotId;
            end
        end
    end

    isCsoSlot = ~isinf(k_xi) && mod(slotId,k_xi)==0;
    if ~isCsoSlot && ~isempty(readyIds)
        pktId = readyIds(1);
        readyIds(1) = [];
        accumulatedWait(pktId) = accumulatedWait(pktId)+slotId-readySince(pktId);
        if TxCount(pktId) == 0
            firstTxSlot(pktId) = slotId;
        elseif isReTxPrioritised
            numReadyRetries = numReadyRetries-1;
        end
        TxCount(pktId) = TxCount(pktId)+1;

        dueBoundary = slotId+1+decodingDelay;
        decodeIndex = mod(dueBoundary-1,decodeRingLength)+1;
        if decodePacket(decodeIndex) ~= 0
            error('CSO:DecodingEventCollision', 'Two packets were scheduled for the same decoding entry.');
        end
        decodePacket(decodeIndex) = pktId;
        decodeSuccess(decodeIndex) = rand < 1-p_ir(TxCount(pktId));
        if recordDetails
            txSlotLog(end+1) = slotId; %#ok<AGROW>
            txPacketLog(end+1) = pktId; %#ok<AGROW>
            txAttemptLog(end+1) = TxCount(pktId); %#ok<AGROW>
            txReadySinceLog(end+1) = readySince(pktId); %#ok<AGROW>
        end
    end

    if numPacketsInSystem < 0 || numPacketsInSystem > maxQueueSize
        error('CSO:PopulationInvariant', 'The maintained population is outside its allowed range.');
    end
    if slotId <= MaxNumSlots
        queueSize(slotId) = numPacketsInSystem;
    end
    slotId = slotId+1;
end

arrSlot = arrSlot(1:numAdmitted);
TxCount = TxCount(1:numAdmitted);
firstTxSlot = firstTxSlot(1:numAdmitted);
completionBoundary = completionBoundary(1:numAdmitted);
discarded = discarded(1:numAdmitted);
accumulatedWait = accumulatedWait(1:numAdmitted);
if any(completionBoundary==0) || any(firstTxSlot==0)
    error('CSO:UnfinishedPackets', 'A packet has no terminal decoding outcome after draining.');
end
totalLatencyInSlotsVector = completionBoundary-arrSlot;
totalLatencyInSlotsVector(discarded) = Inf;
waitLatencyInslotsVector = accumulatedWait;

if recordDetails
    details = struct('generatedPackets',numGenerated,'admittedPackets',numAdmitted, ...
        'overflowDrops',NumDroppedPkts,'successfulPackets',sum(~discarded), ...
        'discardedPackets',sum(discarded),'unfinishedPackets',0, ...
        'drainSlots',max(0,slotId-1-MaxNumSlots),'finalPopulation',numPacketsInSystem, ...
        'decodingDelaySlots',decodingDelay,'feedbackDelaySlots',feedbackDelay, ...
        'arrivalSlot',arrSlot,'firstTxSlot',firstTxSlot,'completionBoundary',completionBoundary, ...
        'populationBeforeArrival',populationBeforeArrival, ...
        'accumulatedWaitingDelay',waitLatencyInslotsVector,'txSlot',txSlotLog, ...
        'txPacketId',txPacketLog,'txAttempt',txAttemptLog,'txReadySince',txReadySinceLog);
end
end
