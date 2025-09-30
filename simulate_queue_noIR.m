function [queueSize,TxCount,totalLatencyInSlotsVector,waitLatencyInslotsVector] = ...
    simulate_queue_noIR...
    (isReTxPrioritised,p_arr,per,feedbackDelay,...
    MaxNumSlots,maxQueueSize,MaxTxPerPkt,isDetArr)
%simulate_queue_noIR Fn. to simulate queue w or w\o feedback delay
%   Basically a markov chain. State observations are done at slot boundary.
%   Arrivals just after and departures just before the slot boundary.
%   Arrivals are processed immediately in the same slot if possible.
%   Failed departures trigger retransmission where the schedule happens in
%   between the departure and the immediate slot boundary. If the feedback
%   delay \delta=0, this schedule is for the immediate next slot; otherwise
%   for the slot \delta after the next slot and the packet is available
%   just before the start of this particular slot. \delta is in slots.


%% Simulation params
% MaxNumSlots = 1e6;
% maxQueueSize = Inf;
% MaxTxPerPkt = Inf;
MaxNumPkts = MaxNumSlots;

%% Initialize
f = p_arr;
c = 1/f;
p_succ = 1-per;
delay_fdbk = feedbackDelay;

%% Consistency Check
if(isDetArr)
    if(floor(c)~=c)
        error("Configured deterministic arrivals but cycle time is not an integer.")
    end
end
if (p_arr>1-per)
    error("Unstable Queue. Change Params")
end


% Preallocate structures
queueSize = zeros(1,MaxNumSlots);
arrSlot = zeros(1,MaxNumPkts);  % Arrived
TxCount = zeros(1,MaxNumPkts);  % Not Transmitted
servSlot = -3 * ones(1,MaxNumPkts);
% >=1: Successfull slotId,
% -1 : Unserved,
% -2 : Discarded.
% -3 : pkt Not arrived
firstTxSlot = zeros(1,MaxNumPkts);

% Other statistics
NumDroppedPkts = 0;
NumIdleSlots = 0;

% Dynamic Params
pktArrId = 0;
pktTxId = 1;
currQueueSize = 0;

txQueue_pktId =[];
reTxPktIdList = -1*ones(1,MaxNumSlots+1);


% Simulate
for slotId = 1:MaxNumSlots
    %% Arrivals
    if(isDetArr)
        isArrival = (mod(slotId,c)==0);
    else
        isArrival = (rand<p_arr);
    end

    if(isArrival)
        if(currQueueSize == maxQueueSize)
            NumDroppedPkts = NumDroppedPkts+1;
        else
            pktArrId = pktArrId + 1;
            arrSlot(pktArrId) = slotId; % Arrived
            firstTxSlot(pktArrId) = -1;
            TxCount(pktArrId) = 0;       % Not Transmitted
            servSlot(pktArrId) = -1;     % -1: Unserved, -2: Discarded

            currQueueSize = currQueueSize + 1;
            txQueue_pktId(end+1) = pktArrId;
        end
    end



    %% Schedule Retransmissions at the start of slot after arrivals
    % % This is commented out/Removed.
    % % Because this model would be inconsistent when feedback delay > 0
    % % When feedback delay is nonzero, the packets in feedback are not
    % % counted towards the queue length and thus in this model, the ReTx
    % % packets will not be conted towards the state observation.

    % % Current model schedules ReTx at the end of slot for the next slot.

    % % If this model is necessary later copy the ReTx section here and
    % % make the following changes
    % % 1) Change reTxPktIdList(slotId + 1) to reTxPktIdList(slotId)
    % % 2) Change "currQueueSize = currQueueSize + 1;" within an
    % % %  "if(delay_fdbk > 0)",
    % % 3) Change "currQueueSize = currQueueSize - 1;" within an
    % % %  "if(delay_fdbk > 0)"
    % % %  in "Transmissions" section under failed transmissions


    %% Transmissions
    if (currQueueSize == 0)
        NumIdleSlots = NumIdleSlots + 1;
    else
        pktTxId = txQueue_pktId(1);
        txQueue_pktId(1)= [];
        if(TxCount(pktTxId) == 0)
            firstTxSlot(pktTxId) = slotId;
        end
        TxCount(pktTxId) = TxCount(pktTxId) + 1;


        if(rand<p_succ) % Successful Transmissions
            servSlot(pktTxId) = slotId;
            currQueueSize = currQueueSize - 1;

        else
            if(TxCount(pktTxId) == MaxTxPerPkt)
                servSlot(pktTxId) = -2; % Discarded
                currQueueSize = currQueueSize - 1;
            else % Schedule for Retransmission
                currQueueSize = currQueueSize - 1;
                % Note for the alternate model where ReTx is avvailable
                % at the start of the slot after/along with Arrivals:
                % In this case the queue size needs to be reduced ONLY
                % when the feedback delay is nonzero. The packet is
                % assumed to be removed from the queue and then later
                % re-added when the necessary retransmission slot
                % arrives. That means, the packets that are currently
                % being fed back are not part of the queue. When feedback delay
                % is zero, the packets are assumed to be not left the
                % queue during a failure and hence no need to reduce
                % the quesize. Read the comment on above section.

                reTxPktIdList(slotId + delay_fdbk + 1) = pktTxId;

                % %                 Debug
                %                 debug_isReTxScheduled(slotId)=1;
            end
        end
    end



    %% Schedule Retransmissions at the end of the slot
    % % Retransmissions are scheduled at the end of the slot.
    % % For eg. with feedbackdelay=0 ReTx at slot n is scheduled at the
    % % very end of the slot n-1 after the transmission is completed and if
    % % it is failed. Note that queuelength is measured @ slot boundary.

    % By design there might be collissions between new packets and ReTx
    % ReTx is prioritized
    % By design there wont be collissions between two retransmissions.

    if(reTxPktIdList(slotId+1) ~= -1)
        currQueueSize = currQueueSize + 1;
        % When the feedback delay is nonzero, the packet is
        % assumed to be removed from the queue and then later
        % re-added here at the retransmission slot.
        if(isReTxPrioritised == true)
            temp = [reTxPktIdList(slotId+1) txQueue_pktId];
            txQueue_pktId = temp;
        else
            txQueue_pktId(end+1) = reTxPktIdList(slotId+1);
        end
    end


    %% Observe the state by measuring the queue length
    queueSize(slotId) = currQueueSize;

    % % Debug
    %     if(slotId>1)
    %     x = (slot(slotId-1).queueSize)+1;
    %     y = (slot(slotId).queueSize)+1;
    %     debug_TransProbMat(x,y) =  debug_TransProbMat(x,y)+1;
    %     end
end

% Queue size is computed at the slot boundary. This is just after the
% departures and just before the arrivals.

isPktUnserved = (servSlot == -1);
isPktNotArrived = (servSlot == -3);
isPktInvalid = (isPktUnserved + isPktNotArrived)>0;
servSlot(isPktInvalid) = [];
TxCount(isPktInvalid) = [];
arrSlot(isPktInvalid) = [];
firstTxSlot(isPktInvalid) = [];

isPktDiscarded = (servSlot == -2);
servSlot(isPktDiscarded) = Inf;

totalLatencyInSlotsVector = 1+ (servSlot - arrSlot);
waitLatencyInslotsVector = (firstTxSlot - arrSlot);

%  % Debug
% pi0 = sum(queueSizeVector==0)/length(queueSizeVector);
% f=p_arr;
% p=per;
% p_noTx =pi0*(1-f);
% [p p*(1-p_noTx) sum(debug_isReTxScheduled)/length(debug_isReTxScheduled)]

end


