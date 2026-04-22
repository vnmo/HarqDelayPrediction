function [probDistr] = getWaitPmf(pVector,f,M,Qmax)
%GETWAITPMF Function to get distribution of wait delay with IR
%   Function takes in probability of error vector of 4 transmissions and
%   the arrival frequency to return the wait distribution. Conditional Wait
%   delay probability at delay = k for a given quelength q is obtained by 
%   the function getWaitProbDensity. Quelength probabilities are obtained
%   by the function getQLenDistribution. Using these two the conditional
%   distribution is marginalised to get the wait distribution.

qLenPmf = getQLenPmf(pVector,f,M,Qmax);
% qLenPmf(qLenPmf<eps)=[];
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