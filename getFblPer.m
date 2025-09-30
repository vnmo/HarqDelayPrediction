function [expectedBler] = getFblPer(iterLen,snrRcvd_avg,specEff,blockLen,V,hvar)

blerArray = zeros(1,iterLen);
fadingCoeff = exprnd(hvar,1,iterLen);

for iter = 1:iterLen
    snrRcvd_inst = (snrRcvd_avg/hvar)*fadingCoeff(iter);
    % snrRcvd_inst = snrTX_avg*fadingCoeff(iter);

    cap_perRe = log2(1+snrRcvd_inst);
    blerArray(iter) = qfunc((cap_perRe - specEff)*sqrt(blockLen/V));
end


expectedBler = mean(blerArray);

end