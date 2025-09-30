function [prob] = getWaitProbability(pdfPoint,queLength,pVector,M)
%GETWAITPROBABILITY Fn. to compute the probability of q packets in k slots.
% Computed using multiple functions for IR and non IR case.
% Non-IR case is not used anymore but used by IR with equal probabilities
% Function also contains unused legacy codes which are deprecated.

k = pdfPoint; %in slots
q = queLength;


if( k<q || k>M*q )
    prob = 0;
    return
elseif length(pVector) == M
    % prob = getprob_IR_deprecated(k,q,pVector); % DEPRECATED
    prob = getprob_IR(k,q,pVector,M,M);
    return
elseif length(pVector) == 1
    % prob = getprob_noIR(k,q,pVector); % DEPRECATED
    pVector = ones(1,M)*pVector;
    prob = getprob_IR(k,q,pVector,M,M);
    return
else
    error("incompatible Probability Vector for IR");
end

end



function [prob] = getprob_IR(k,q,pVector,M,MaxTx_IR)
%getprob_IRfinal Fn. to compute the probability of q packets in k slots.
% % Start with max allowed Tx M, and max possible packets with M Tx denoted
% by N. We will compute probability of all these packets recursively by
% reducing first n=N:-1:0 in the outer loop and within each loop reducing
% the max allowed Tx m recursively from m = Max possible M. (=4 typically).
% Cases with, for example, no packets with M tx would be taken care of the
% loop with n=0, which points to no packets with M tx.

% MaxTx_IR = M in the original recursion loop. MaxTx_IR stays unchanged
% while M is changed within the recursion. this is because we need to
% handle the discard case of MaxTx_IR, and thus cannot loose the variable.

% % Exit clause
% These conditions act both as an early exit as well as a regular exit
% while using as a recursion.
if(k==q)
    % k=q => all success. Return probability
    prob = (1-pVector(1))^q;
    return
elseif( k<q || k>M*q )
    prob = 0;   % k out of possible range => prob = 0.
    return
end

prob = 0;

% Mth Tx of any packet is special case due to potential discard.

N = min(floor((k-q)/(M-1)),q);  % Max num packets with M attempts
% N can never exceed the total num of packets, hence the minimum with q.

% % Note:
% There are cases where N packets of m attempts in invalid. For example
% k-7,q=3,m=2 => N=3. But 3 packets of 2 slots each adds up to only 6 but
% not 7 and hence only possibilities are 1 and 2 but not 3 in the below
% loop with loop variable n. This however, returns a probability
% of zero and hence 'auto managed';
for n = 0:N
    % Loop for all possible seq with n<=N pkts with m Tx attempts.
    % n is the number of pkts with m-1 fails followed by a succ/discard
    % n = 0 => no packets with M transmission and need to check with M-1.

    % Find number ways n such pkts can be chosen from q packets
    numSequences =nchoosek(q , n);

    % In each of these sequences, there are m-1 failures Find prob.
    seqProb_failures  = (prod(pVector(1:M-1)));

    % In these sequences, there is also the final transmission which is
    % a success if m < M. If m=M, it can also be a discard.
    % Both succ & failure at M^th Tx gives same sequence length.
    % Thus corresponding probability = p(succ)+p(fail)=1.
    if(M==MaxTx_IR) % M of recursion = original M stored in MaxTx_IR
        seqProb_success =1;
    else % m != M
        seqProb_success = (1-pVector(M));
    end
    seqProb = seqProb_failures * seqProb_success; % Prob of 1 such seq.
    seqProb_total = seqProb^n; % Total prob. of n such sequences. iid.

    % With the prob. of all n packets computed, we need to find the
    % probability of all possible subsequences of all the other packets
    % with m-1 tranmsissions or less, recursively.
    %% Recursion to find the probability of the subsequences.
    subSeqProb = getprob_IR(k-M*n,q-n,pVector,M-1,MaxTx_IR);

    % Total probability
    prob = prob + numSequences * seqProb_total * subSeqProb;
end
end


function [prob] = getprob_noIR(k,q,p)
if(k==q)
    prob = (1-p)^q;
    return
end

%% First we consider cases without a 4th transmission.

numSolWithoutTx4 = getSolutions(k, q , 3);
% Total ways q packets of max 3 Tx each can make up a delay of k slots

prob = numSolWithoutTx4 * (1-p)^(q) * p^(k-q);
% Total combined Probability of these transmissions
% There are q packets, so total q success and remaining k-q failures
% Works with this simple expression due to iid. (note: no IR here)

%% Now we deal where there is a 4th Tx, success, or discard.
% Let n be the number of transmissions out of q packets
% that went on to 4th Transmissions.

% count_debug = numSolWithoutTx4;
len = floor((k-q)/3);
for n = 1:len
    k1 = k-4*n;
    q1 = q-n;
    prob = prob + nchoosek(q , n) *  (p^3)^n ...
        * getSolutions(k1, q1 , 3)* (1-p)^(q1) * p^(k1-q1);
    %     count_debug = count_debug + nchoosek(q , n) * getSolutions(k1, q1 , 3) ;

end
% count_debug
end


function [numCombinations] = getSolutions(k,q,maxTx)
% The number of solutions to the equation:
% 𝑥1 + 𝑥2 + ... 𝑥q = k, where 1 ≤ 𝑥i ≤ maxTx for all 1 ≤ i ≤ q
M = maxTx;
if(k==q)
    numCombinations=1;
    return
end
len = floor((k-q)/M);
numCombinations = 0;
for i=0:len
    numCombinations = numCombinations + nchoosek(q , i) * (-1)^i *  nchoosek(k-1-M*i , q-1);
end
end


function [prob] = getprob_IR_deprecated(k,q,pVector)

% M represent the maximum allowed number of transmissions in the system
if(k==q)
    % k=q => all success. Return probability
    prob = (1-pVector(1))^q;
    return
end

prob = 0;
M=4;
% Mth Tx of any packet is special case due to potential discard.
% Start from this case and then recursively find the complete probability.
% Loop over all seq where there is at least one packet with m attempts, and
% all other packets taking less than m attempts to cover all possibilities.
% m=1 (equivalently k=q) is not included as it is already considered.
for m=[4 3 2]
    N = min(floor((k-q)/(m-1)),q);
    % Max num packets with m attempts
    % Can never exceep the total num of packets, hence the minimum.

    % % Note:
    % There are cases where N packets of m attempts in invalid. For
    % example k-7,q=3,m=2 => N=3. But 3 packets of 2 slots each adds up to
    % only 6 but not 7 and hence only possibilities are 1 and 2 but not 3
    % in the below lop of "for n=1:N". This however returns a probability
    % of zero and hence 'auto managed';
    % % % % %

    for n = 1:N
        % Loop for all possible seq with n<=N pkts with m Tx attempts.
        % n is the number of pkts with m-1 fails followed by a succ/discard
        % n =0 TODO Explain

        % Find number ways n such pkts can be chosen from q packets
        numSequences =nchoosek(q , n);

        % In each of these sequences, there are m-1 failures Find prob.
        seqProb_failures  = (prod(pVector(1:m-1)));

        % In these sequences, there is also the final transmission which is
        % a success if m < M. If m=M, it can also be a discard.
        % Both succ & failure at M^th Tx gives same sequence length.
        % Thus corresponding probability = p(succ)+p(fail)=1.
        if(m==M) % m = M
            seqProb_success =1;
        else % m != M
            seqProb_success = (1-pVector(m));
        end
        seqProb = seqProb_failures * seqProb_success; % Prob of 1 such seq.
        seqProb_total = seqProb^n; % Total prob. of n such sequences. iid.

        % With the prob. of all n packets computed, we need to find the
        % probability of all possible subsequences of all the other packets
        % with m-1 tranmsissions or less, recursively.
        subSeqProb = getPdf_IR_subseq_deprecated(k-m*n,q-n,m-1,pVector);

        % Total probability
        prob = prob + numSequences * seqProb_total * subSeqProb;
    end
end
end


function [prob] = getPdf_IR_subseq_deprecated(k,q,m,pVector)
% Function to dind the probability of the subsequence.
% The subseq contains q packets and k slots. Max allowed Tx per pkt is m.
prob = 0;
if(k==q)
    prob = (1-pVector(1))^q;    % k=q => all success. Return probability
    return
elseif( k<q || k>m*q )
    prob = 0;   % k out of possible range => prob = 0.
    return
end

% MaxTx=m;
%  for m=MaxTx:-1:2
% Let N is Max #pkts with m-1 failures and a success.
% No discards possible in any subsequence by construction.
% Discards considered already in the main function.
% m = 1 would be already caught in the (k==q) or (k<q || k>m*q) before.
N = floor((k-q)/(m-1));


% Loop for all possible seq with n<=N pkts with m Tx attempts.
%TODO why not starting from 1
for n = 0:N
    numSequences = nchoosek(q , n);
    seqProb_total = (prod(pVector(1:m-1))*(1-pVector(m)))^n;
    subSeqProb = getPdf_IR_subseq_deprecated(k-m*n,q-n,m-1,pVector);

    prob = prob + numSequences * seqProb_total * subSeqProb;
end
%  end
end

