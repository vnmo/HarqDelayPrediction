function [prob] = getWaitProbability( ...
    pdfPoint, queLength, pVector, M)

k = pdfPoint;
q = queLength;

if k < q || k > M * q
    prob = 0;
    return;
end

if length(pVector) == 1
    pVector = repmat(pVector, 1, M);
elseif length(pVector) >= M
    pVector = pVector(1:M);
else
    error("pVector must contain at least M elements.");
end

pVector = pVector(:).';

prob = getprob_IR(k, q, pVector, M, M);

end


function [prob] = getprob_IR( ...
    k, q, pVector, M, MaxTx_IR)

if k == q
    if MaxTx_IR == 1
        prob = 1;
    else
        prob = (1 - pVector(1))^q;
    end
    return;
elseif k < q || k > M * q
    prob = 0;
    return;
end

prob = 0;

N = min( ...
    floor((k - q) / (M - 1)), ...
    q);

for n = 0:N

    numSequences = nchoosek(q, n);

    seqProb_failures = prod(pVector(1:M-1));

    if M == MaxTx_IR
        seqProb_success = 1;
    else
        seqProb_success = 1 - pVector(M);
    end

    seqProb = ...
        seqProb_failures * seqProb_success;

    seqProb_total = seqProb^n;

    subSeqProb = getprob_IR( ...
        k - M * n, ...
        q - n, ...
        pVector, ...
        M - 1, ...
        MaxTx_IR);

    prob = prob + ...
        numSequences * seqProb_total * subSeqProb;
end

end
