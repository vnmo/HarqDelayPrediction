function [expectedBlerVector] = getFblPerVector_Harq(iterLen,snrRcvd_avg,specEff,blockLen,V,hvar,M_Harq)
error("TO complete. Doublecheck: update required for the comment by Mikael")

blerArray = zeros(iterLen,M_Harq);
fadingCoeff = exprnd(hvar,4,iterLen);

for iter = 1:iterLen
    snrRcvd_inst = (snrRcvd_avg/hvar)*fadingCoeff(:,iter);
    % snrRcvd_inst = snrTX_avg*fadingCoeff(iter);
   
    for m=1:M_Harq
        cap_perRe = mean(log2(1+snrRcvd_inst(1:m)));
        blerArray(iter,m) = qfunc((cap_perRe - (specEff/m))*sqrt(m * blockLen/V));
    end
end

expectedBlerVector = mean(blerArray);

end

% blerArray = zeros(iterLen,M_Harq);
% fadingCoeff = exprnd(hvar,1,iterLen);
% 
% for iter = 1:iterLen
%     snrRcvd_inst = (snrRcvd_avg/hvar)*fadingCoeff(iter);
%     % snrRcvd_inst = snrTX_avg*fadingCoeff(iter);
% 
%     cap_perRe = log2(1+snrRcvd_inst);
%     for m=1:M_Harq
%         blerArray(iter,m) = qfunc((cap_perRe - (specEff/m))*sqrt(m * blockLen/V));
%     end
% end
% 
% expectedBlerVector = mean(blerArray);
% 
% end