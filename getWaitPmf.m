function [probDistr] = getWaitPmf(pVector,f,M,Qmax)

qLenPmf = getQLenPmf(pVector,f,M,Qmax);
p_emptyQ = qLenPmf(1);
qLenPmf(1) = [];
Lq=length(qLenPmf);
Lw = 4*Lq+1;
waitDistribution =zeros(1,Lw);

for qLen=1:Lq
    condWaitDistribution = zeros(1,Lw);
    for k =qLen:Lw
        wdensity = getWaitProbability(k,qLen,pVector,M);
        if (wdensity<eps)
            break;
        end
        condWaitDistribution(k) = wdensity;
    end
    waitDistribution = waitDistribution + qLenPmf(qLen)*condWaitDistribution;
end
probDistr = [p_emptyQ waitDistribution];
end