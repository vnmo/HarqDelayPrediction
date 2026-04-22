function [prob] = getprob_IR(k,q,pVector,M,MaxTx_IR)
    if(k==q)
        prob = (1-pVector(1))^q;
        return
    elseif( k<q || k>M*q )
        prob = 0;
        return
    end

    prob = 0;
    N = min(floor((k-q)/(M-1)),q);

    for n = 0:N
        numSequences = nchoosek(q , n);
        seqProb_failures  = (prod(pVector(1:M-1)));

        if(M==MaxTx_IR)
            seqProb_success = 1;
        else
            seqProb_success = (1-pVector(M));
        end

        seqProb = seqProb_failures * seqProb_success;
        seqProb_total = seqProb^n;

        subSeqProb = getprob_IR(k-M*n, q-n, pVector, M-1, MaxTx_IR);

        prob = prob + numSequences * seqProb_total * subSeqProb;
    end
end