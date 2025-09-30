function [QLenPmf] = getQLenPmf(pVector,f,M,Qmax)
%GETQLENPMF Function to get steady state Q probabilities with IR.
%   Function takes in probability of error vector of 4 transmissions and
%   the arrival frequency to return the steady state queue probabilities.
%   For this an assumption of max HARQ processes of 16 is made. Markov
%   transition matrix is constructed with each state corresponding to a
%   (quelength, transmission number) tuple, which is made into a continuous
%   number from 0 corresponding to (0,0) and 4*(16+1) corresponding to
%   (16,4). Steady state probability is the eigen vector corresponding to
%   unit eigenvalue of the transpose of the transition probability matrix.

if(length(pVector) ==1)
    p(1:4)=pVector;
elseif(length(pVector) ==4)
    p=pVector;
else
    error("incompatible");
end


%% Deprecated code segment. 
% % A better implementation in sync with the writing and easier to interpret
% % is avalaible below, with potential cost of small runtime increase.

% L = (16+1)*4; %(Qmax+1)M
% P = zeros(L,L);
% 
% P(1,1)=1-f*p(1);
% P(1,6)=f*p(1);
% P(2,2)=1;
% P(3,3)=1;
% P(4,4)=1;
% 
% % i corresponds to queueSize. Max QueueSize is 16, and minimum is 0
% for i=1:16
%     % j corresponds to Total Tx attempts completed by next packet 
%     % j=0 => Fresh packet
%     % j=1 => Second Tx
%     % j=2 => Third Tx
%     % j=3 => Fourth Tx
%     for j=0:3
%         r = i*4+(j+1);
%         if(j==3)     
%             P(r,r-4-j)  = 1-f;
%             P(r,r-j)    = f;
%         elseif(i==16)
%             P(r,r+1)    = p(j+1);
%             P(r,r-4-j)  = (1-f)*(1-p(j+1));
%             P(r,r-j)    = f*(1-p(j+1));
%         else
%             P(r,r+1)    = (1-f)*p(j+1);
%             P(r,r+5)    = f*p(j+1);
%             P(r,r-4-j)  = (1-f)*(1-p(j+1));
%             P(r,r-j)    = f*(1-p(j+1));
%         end
%     end
% end
% 
% P1=P;

L = (Qmax+1)*M;
p_prime = 1 - p;  % complement of p
f_prime = 1 - f;  % complement of f

P = zeros(L,L);
for i = 2:M
    P(i,i) = 1; % (0,m), m>=2l unreachable states.
end

% Transition (0,1) -> (0,1)
q = 0; m = 1; s1 = q * M + m;
P(s1, s1) = 1 - f * p(1);

% Transition (0,1) -> (1,2)
q = 1; m = 2; s2 = q * M + m;
P(s1, s2) = f * p(1);

% Iterate through all values of q and m
for q = 1:Qmax
    for m = 1:M
        s1 = q * M + m;  % Current state (q,m)

        if q < Qmax && m < M
            % Next state (q,m+1)
            s2 = q * M + (m + 1);  
            P(s1, s2) = f_prime * p(m);  

            % Next state (q+1,m+1)
            s2 = (q + 1) * M + (m + 1);
            P(s1, s2) = f * p(m);  
        end

        if m < M
            % Next state (q,1)
            s2 = q * M + 1;   
            P(s1, s2) = f * p_prime(m);  

            % Next state (q-1,1)
            s2 = (q - 1) * M + 1;   
            P(s1, s2) = f_prime * p_prime(m);  
        end

        if m == M
            % Next state (q,1)
            s2 = q * M + 1;
            P(s1, s2) = f;  

            % Next state (q-1,1)
            s2 = (q - 1) * M + 1; 
            P(s1, s2) = f_prime;  
        end

        if q == Qmax && m < M
            % Next state (Qmax,m+1)
            s2 = q * M + (m + 1);
            P(s1, s2) = p(m); 
        end
    end
end


[ V , D ] = eig( P' ); % Find eigenvector/eigenvalue

for i =1:size(D,1)
    % Double confirming and choosing  the eigenvector of unit eigenvalue
    if (D(i,i)>1-sqrt(eps) && D(i,i)<1+sqrt(eps)) 
        vec = V(:,i); 
        break;
    end
end

% Normalise to 1 to represent probability
vec = vec/sum(vec);

% Marginalise over m for each q to get probability of queuelengths
QLenPmf=sum(reshape(vec,M,Qmax+1));

end