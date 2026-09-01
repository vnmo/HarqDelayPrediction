function [queueSize,TxCount,totalLatencyInSlotsVector,waitLatencyInslotsVector] = ...
    simulate_queue_noIR...
    (isReTxPrioritised,p_arr,per,feedbackDelay,...
    MaxNumSlots,maxQueueSize,MaxTxPerPkt,isDetArr)


MaxNumPkts = MaxNumSlots;

f = p_arr;
c = 1/f;
p_succ = 1-per;
delay_fdbk = feedbackDelay;

if(isDetArr)
    if(floor(c)~=c)
        error("Configured deterministic arrivals but cycle time is not an integer.")
    end
end
if (p_arr>1-per)
    error("Unstable Queue. Change Params")
end


queueSize = zeros(1,MaxNumSlots);
arrSlot = zeros(1,MaxNumPkts);  % Arrived
TxCount = zeros(1,MaxNumPkts);  % Not Transmitted
servSlot = -3 * ones(1,MaxNumPkts);
firstTxSlot = zeros(1,MaxNumPkts);

NumDroppedPkts = 0;
NumIdleSlots = 0;

pktArrId = 0;
pktTxId = 1;
currQueueSize = 0;

txQueue_pktId =[];
reTxPktIdList = -1*ones(1,MaxNumSlots+1);


for slotId = 1:MaxNumSlots
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







    if (currQueueSize == 0)
        NumIdleSlots = NumIdleSlots + 1;
    else
        pktTxId = txQueue_pktId(1);
        txQueue_pktId(1)= [];
        if(TxCount(pktTxId) == 0)
            firstTxSlot(pktTxId) = slotId;
        end
        TxCount(pktTxId) = TxCount(pktTxId) + 1;


        if(rand<p_succ)
            servSlot(pktTxId) = slotId;
            currQueueSize = currQueueSize - 1;

        else
            if(TxCount(pktTxId) == MaxTxPerPkt)
                servSlot(pktTxId) = -2; % Discarded
                currQueueSize = currQueueSize - 1;
            else
                currQueueSize = currQueueSize - 1;

                reTxPktIdList(slotId + delay_fdbk + 1) = pktTxId;

            end
        end
    end





    if(reTxPktIdList(slotId+1) ~= -1)
        currQueueSize = currQueueSize + 1;
        if(isReTxPrioritised == true)
            temp = [reTxPktIdList(slotId+1) txQueue_pktId];
            txQueue_pktId = temp;
        else
            txQueue_pktId(end+1) = reTxPktIdList(slotId+1);
        end
    end


    queueSize(slotId) = currQueueSize;

end


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


end

