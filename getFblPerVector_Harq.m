
function [expectedBlerVector] = getFblPerVector_Harq(iterLen,snrRcvd_avg,specEff,blockLen,V,hvar,M_Harq)

blerArray = zeros(iterLen,M_Harq);
fadingCoeff = exprnd(hvar,4,iterLen);

for iter = 1:iterLen
    snrRcvd_inst = (snrRcvd_avg/hvar)*fadingCoeff(:,iter);
   
    for m=1:M_Harq
        cap_perRe = mean(log2(1+snrRcvd_inst(1:m)));
        blerArray(iter,m) = qfunc((cap_perRe - (specEff/m))*sqrt(m * blockLen/V));
    end
end

expectedBlerVector = mean(blerArray);

end

