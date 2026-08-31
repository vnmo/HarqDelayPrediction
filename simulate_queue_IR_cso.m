function [queueSize, TxCount, totalLatencyInSlotsVector, waitLatencyInslotsVector, NumDroppedPkts] = ...
    simulate_queue_IR_cso( ...
    isReTxPrioritised, p_arrival, perVector, feedbackDelay, k_xi, ...
    MaxNumSlots, maxQueueSize, MaxTxPerPkt, isDetArr)
%SIMULATE_QUEUE_IR_CSO Slot-by-slot IR-HARQ queueing simulation with CSO.
% This function samples arrivals and transmission outcomes. It does not
% construct or call the Proposed model's (q,m) transition matrix.
% k_xi is the CSP length in slots, not the current phase. The current phase
% is mod(slotId,k_xi); phase zero blocks data transmission for one slot.
% Decoding delay is fixed to one slot. A failed packet becomes eligible
% again after 1+decodingDelay+feedbackDelay slots from its transmission.
% Other ready packets may transmit while that packet waits for feedback.
% Packet lists and scheduled retry slots store the simulation state;
% no finite HARQ-process pool is allocated here.
%
% Per-slot order: arrival/admission, retry availability, CSO check, then
% at most one data transmission. The branches below describe each update.
% totalLatencyInSlotsVector ends at decoding completion, not final ACK.
% waitLatencyInslotsVector measures only arrival-to-first-transmission
% waiting; it is not the auxiliary waiting metric returned by getWaitPmfCSO.

%% Simulation params
MaxNumPkts = MaxNumSlots;

p_ir = perVector;
if length(p_ir) ~= MaxTxPerPkt
    error("Size of probability vector should be " + MaxTxPerPkt)
end

expTxPerPkt = 0;
for j = 1:MaxTxPerPkt
    expTxPerPkt = expTxPerPkt + j * (1 - p_ir(j)) * prod(p_ir(1:j-1));
end
expTxPerPkt = expTxPerPkt + MaxTxPerPkt * prod(p_ir);

%% Initialize
f = p_arrival;
delay_fdbk = feedbackDelay;
decodingDelay = 1;               % fixed decoding delay (one slot)
c = 1 / f;

%% Consistency check
if isDetArr
    if floor(c) ~= c
        error("Configured deterministic arrivals but cycle time is not an integer.")
    end
end

% Effective service-rate check with CSO
effServiceRate = (1 - 1 / k_xi) / expTxPerPkt;
if p_arrival > effServiceRate
    error("Unstable Queue. Change Params")
end

%% Preallocate structures
% Per-packet state is stored by packet ID, rather than by a DTMC index.
% servSlot=-1 covers both ready packets and packets awaiting a retry.
% A positive servSlot stores scheduled decoding completion after success.
queueSize   = zeros(1, MaxNumSlots);
arrSlot     = zeros(1, MaxNumPkts);      % Arrival slot
TxCount     = zeros(1, MaxNumPkts);      % Number of transmissions
servSlot    = -3 * ones(1, MaxNumPkts);  % >=1: completion slot
                                         % -1 : unserved
                                         % -2 : discarded
                                         % -3 : not arrived
firstTxSlot = zeros(1, MaxNumPkts);      % First transmission slot

%% Statistics
NumDroppedPkts = 0;
NumIdleSlots   = 0;
NumBlockedSlots = 0; %#ok<NASGU>

%% Dynamic variables
% txQueue_pktId holds ready packets; currQueueSize is its length.
% numPacketsInSystem also counts packets awaiting a retry. The counter
% decreases when success/discard is sampled, although successful latency
% is recorded at the later decoding-completion slot in servSlot.
pktArrId = 0;
currQueueSize = 0;
numPacketsInSystem = 0;
txQueue_pktId = [];

% reTxPktIdList(t) = packet that becomes available at beginning of slot t
reTxPktIdList = -1 * ones(1, MaxNumSlots + decodingDelay + delay_fdbk + 2);

%% Simulation loop
for slotId = 1:MaxNumSlots

    %% Arrivals
    if isDetArr
        isArrival = (mod(slotId, c) == 0);
    else
        isArrival = (rand < p_arrival);
    end

    if isArrival
        if numPacketsInSystem >= maxQueueSize
            % Full-buffer branch: reject the arrival; stored packets remain.
            NumDroppedPkts = NumDroppedPkts + 1;
        else
            % Admission branch: new packet is ready at the queue tail,
            % with TxCount=0 and servSlot=-1 until success or discard.
            pktArrId = pktArrId + 1;
            arrSlot(pktArrId) = slotId;
            firstTxSlot(pktArrId) = -1;
            TxCount(pktArrId) = 0;
            servSlot(pktArrId) = -1;

            currQueueSize = currQueueSize + 1;
            numPacketsInSystem = numPacketsInSystem + 1;
            txQueue_pktId(end+1) = pktArrId;
        end
    end

    %% Retransmissions become available at the start of this slot
    if reTxPktIdList(slotId) ~= -1
        % Awaiting retry -> ready. The packet is already admitted, so
        % numPacketsInSystem does not increase again.
        currQueueSize = currQueueSize + 1;

        if isReTxPrioritised
            % A newly eligible retry is inserted at the queue head.
            txQueue_pktId = [reTxPktIdList(slotId), txQueue_pktId];
        else
            % Without retry priority, insert it at the queue tail.
            txQueue_pktId(end+1) = reTxPktIdList(slotId);
        end
    end

    %% CSO check: one blocked slot every k_xi slots
    % The clock supplies the CSP phase; no separate phase array is needed.
    % mod(slotId,k_xi)=0 -> blocked; phases 1,...,k_xi-1 -> data slot.
    isCsoSlot = (mod(slotId, k_xi) == 0);

    %% Transmission
    if isCsoSlot
        % Blocked branch: no packet is transmitted and no TxCount changes.
        % Arrivals and retry availability have already been processed.
        NumBlockedSlots = NumBlockedSlots + 1;

    elseif currQueueSize == 0
        % Idle branch: no packet is ready; future retries may still exist.
        NumIdleSlots = NumIdleSlots + 1; %#ok<NASGU>

    else
        % Data branch: transmit the head packet and increment its count.
        pktTxId = txQueue_pktId(1);
        txQueue_pktId(1) = [];

        if TxCount(pktTxId) == 0
            firstTxSlot(pktTxId) = slotId;
        end

        TxCount(pktTxId) = TxCount(pktTxId) + 1;

        pSuccess = 1 - p_ir(TxCount(pktTxId));
        if rand < pSuccess
            % Success branch: ready -> scheduled decoding completion.
            % Release the packet count now; add decoding to its latency.
            servSlot(pktTxId) = slotId + decodingDelay;
            currQueueSize = currQueueSize - 1;
            numPacketsInSystem = numPacketsInSystem - 1;

        else
            if TxCount(pktTxId) == MaxTxPerPkt
                % Final failure branch: ready -> discarded (servSlot=-2).
                servSlot(pktTxId) = -2;
                currQueueSize = currQueueSize - 1;
                numPacketsInSystem = numPacketsInSystem - 1;
            else
                % Retry branch: ready -> awaiting retry (servSlot stays -1).
                % It leaves the ready queue but remains admitted.
                % futureSlot includes transmission, decoding, and feedback.
                currQueueSize = currQueueSize - 1;
                futureSlot = slotId + decodingDelay + delay_fdbk + 1;
                if futureSlot <= length(reTxPktIdList)
                    reTxPktIdList(futureSlot) = pktTxId;
                end
            end
        end
    end

    %% Observe queue size at slot boundary
    queueSize(slotId) = numPacketsInSystem;
end

%% Post-processing
% Omit unresolved packets and unused entries from the returned latency
% arrays. Admission drops are returned separately in NumDroppedPkts.
isPktUnserved   = (servSlot == -1);
isPktNotArrived = (servSlot == -3);
isPktInvalid    = (isPktUnserved + isPktNotArrived) > 0;

servSlot(isPktInvalid) = [];
TxCount(isPktInvalid) = [];
arrSlot(isPktInvalid) = [];
firstTxSlot(isPktInvalid) = [];

isPktDiscarded = (servSlot == -2);
% Discards have infinite latency and therefore violate every finite target.
servSlot(isPktDiscarded) = Inf;

% Include the transmission slot and successful decoding, but no final ACK.
totalLatencyInSlotsVector = 1 + (servSlot - arrSlot);
% This diagnostic output excludes waiting between retransmissions.
waitLatencyInslotsVector  = firstTxSlot - arrSlot;

end
