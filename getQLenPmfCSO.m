% Deterministic
% function [QLenPmf] = getQLenPmfCSO(pVector,f,M,Qmax,CSO)
% Deterministic periodic ON/OFF schedule Markov chain with state (q,m,k).
% q=0..Qmax, m=1..M, k=0..T-1, k'=(k+1) mod T.
% %
% State index (MATLAB 1-index):
%   s = k*(Qmax+1)*M + q*M + m
% %
% CSO:
%   - If vector: CSO(k)=xi(k) in {0,1} (1=ON, 0=OFF), length T
%   - If scalar: interpret CSO = 1/T_off (e.g., 1/80) => build T=round(1/CSO),
%                xi=1 except one OFF slot xi=0
% 
%     ---- p vector ----
%     if length(pVector)==1
%         p = repmat(pVector,1,M);
%     elseif length(pVector)==M
%         p = pVector(:).';
%     else
%         error("incompatible pVector length");
%     end
%     p_prime = 1 - p;
%     f_prime = 1 - f;
% 
%     ---- schedule xi(k) ----
%     if isempty(CSO)
%         error("CSO must be provided.");
%     end
% 
%     if isscalar(CSO)
%         T = round(1/CSO);
%         if T < 1
%             error("Scalar CSO is invalid.");
%         end
%         xiVec = ones(1,T);   ON
%         xiVec(1) = 0;        one OFF per period (원하면 위치 바꿔도 됨)
%     else
%         xiVec = CSO(:).';
%         T = length(xiVec);
%         if any((xiVec~=0) & (xiVec~=1))
%             error("CSO schedule vector must contain only 0/1.");
%         end
%     end
% 
%     Lq = Qmax+1;
%     block = Lq*M;            states per k
%     L = block*T;             total states
% 
%     index helper (q in 0..Qmax, m in 1..M, k0 in 0..T-1)
%     idx = @(q,m,k0) (k0*block + q*M + m);
% 
%     sparse triplets
%     I = []; J = []; V = [];
%     function add(i,j,v)
%         if v~=0
%             I(end+1,1)=i;
%             J(end+1,1)=j;
%             V(end+1,1)=v;
%         end
%     end
% 
%     ---- Build P row by row ----
%     for k0 = 0:T-1
%         xi = xiVec(k0+1);
%         kp = mod(k0+1,T);
% 
%         unreachable (0,m>=2,k): force (0,m,k)->(0,m,k') prob 1
%         for m = 2:M
%             s = idx(0,m,k0);
%             add(s, idx(0,m,kp), 1);
%         end
% 
%         ---------- (0,1,k) ----------
%         s = idx(0,1,k0);
%         rowSum = 0;
% 
%         (0,1,k)->(0,1,k')
%         v = 1 - (1-xi)*f - xi*f*p(1);
%         add(s, idx(0,1,kp), v); rowSum = rowSum + v;
% 
%         (0,1,k)->(1,1,k')
%         if Qmax >= 1
%             v = (1-xi)*f;
%             add(s, idx(1,1,kp), v); rowSum = rowSum + v;
% 
%             (0,1,k)->(1,2,k')
%             if M >= 2
%                 v = xi*f*p(1);
%                 add(s, idx(1,2,kp), v); rowSum = rowSum + v;
%             end
%         end
% 
%         just in case of rounding
%         if abs(rowSum-1) > 1e-12
%             add(s, idx(0,1,kp), 1-rowSum);
%         end
% 
%         ---------- q = 1..Qmax ----------
%         for q = 1:Qmax
%             for m = 1:M
%                 s1 = idx(q,m,k0);
%                 rowSum = 0;
% 
%                 ----- boundary q==Qmax: use your bottom block -----
%                 if q == Qmax
%                     if xi == 0
%                         (Qmax,m,k)->(Qmax,m,k') with prob 1
%                         add(s1, idx(Qmax,m,kp), 1);
%                         rowSum = 1;
%                     else
%                         if m < M
%                             (Qmax,m,k)->(Qmax,m+1,k') prob p_m
%                             v = p(m);
%                             add(s1, idx(Qmax,m+1,kp), v); rowSum=rowSum+v;
% 
%                             (Qmax,m,k)->(Qmax-1,1,k') prob p_m'
%                             v = p_prime(m);
%                             add(s1, idx(Qmax-1,1,kp), v); rowSum=rowSum+v;
%                         else
%                             (Qmax,M,k)->(Qmax-1,1,k') prob 1
%                             add(s1, idx(Qmax-1,1,kp), 1);
%                             rowSum = 1;
%                         end
%                     end
% 
%                     (안전) 행합이 1이 아니면 (Qmax,m,k')-> stay에 넣기
%                     if abs(rowSum-1) > 1e-12
%                         add(s1, idx(Qmax,m,kp), 1-rowSum);
%                     end
%                     continue;
%                 end
% 
%                 ----- interior q < Qmax -----
% 
%                 (q,1,k)->(q,1,k') : (1-xi)f' + xi*f*p1'
%                 if m == 1
%                     v = (1-xi)*f_prime + xi*f*(1-p(1));
%                     add(s1, idx(q,1,kp), v); rowSum=rowSum+v;
%                 end
% 
%                 (q,m,k)->(q,m,k') for m=2..M-1 : (1-xi)f'
%                 if (m >= 2) && (m <= M-1)
%                     v = (1-xi)*f_prime;
%                     add(s1, idx(q,m,kp), v); rowSum=rowSum+v;
%                 end
% 
%                 (q,m,k)->(q+1,m,k') for m=1..M-1 : (1-xi)f
%                 if m <= M-1
%                     v = (1-xi)*f;
%                     add(s1, idx(q+1,m,kp), v); rowSum=rowSum+v;
%                 end
% 
%                 ON part for m=1..M-1
%                 if m <= M-1
%                     (q,m+1) : xi*f'*p_m
%                     v = xi*f_prime*p(m);
%                     add(s1, idx(q,m+1,kp), v); rowSum=rowSum+v;
% 
%                     (q+1,m+1) : xi*f*p_m
%                     v = xi*f*p(m);
%                     add(s1, idx(q+1,m+1,kp), v); rowSum=rowSum+v;
% 
%                     (q,1) : xi*f*p_m'  (주의: 이 행은 m=2..M-1에 해당하는 걸로 표에 적었지만,
%                     m=1일 때는 같은 상태로 가는 항이라 위 self 항에 이미 들어가 있어야 해서 중복 방지)
%                     if m >= 2
%                         v = xi*f*p_prime(m);
%                         add(s1, idx(q,1,kp), v); rowSum=rowSum+v;
%                     end
% 
%                     (q-1,1) : xi*f'*p_m'
%                     v = xi*f_prime*p_prime(m);
%                     add(s1, idx(q-1,1,kp), v); rowSum=rowSum+v;
%                 end
% 
%                 m==M extra rows (your table)
%                 if m == M
%                     v = xi*f;        add(s1, idx(q,1,kp),   v); rowSum=rowSum+v;
%                     v = xi*f_prime;  add(s1, idx(q-1,1,kp), v); rowSum=rowSum+v;
%                     v = (1-xi)*f_prime; add(s1, idx(q,M,kp),   v); rowSum=rowSum+v;
%                     v = (1-xi)*f;       add(s1, idx(q+1,M,kp), v); rowSum=rowSum+v;
%                 end
% 
%                 ---- fill remaining mass to "stay" (q,m,k') ----
%                 if abs(rowSum-1) > 1e-12
%                     add(s1, idx(q,m,kp), 1-rowSum);
%                 end
%             end
%         end
%     end
% 
%     build sparse P (duplicates automatically summed)
%     P = sparse(I,J,V,L,L);
% 
%     ---- steady state ----
%     power iteration is robust for this size
%     pi = ones(L,1)/L;
%     for it = 1:50000
%         pi2 = P.' * pi;
%         pi2 = pi2 / sum(pi2);
%         if norm(pi2-pi,1) < 1e-13
%             pi = pi2;
%             break;
%         end
%         pi = pi2;
%     end
% 
%     ---- marginalize over m,k -> QLenPmf(q) ----
%     pi_mqk = reshape(pi, [M, Lq, T]);   (m, q, k)
%     QLenPmf = squeeze(sum(sum(pi_mqk, 1), 3)); sum over m,k
%     QLenPmf = QLenPmf(:).';
%     QLenPmf
% 
% end


function [QLenPmf] = getQLenPmfCSO(pVector,f,M,Qmax,CSO)
  %GETQLENPMF Function to get steady state Q probabilities with IR.
  %Function takes in probability of error vector of 4 transmissions and
  %the arrival frequency to return the steady state queue probabilities.
  %For this an assumption of max HARQ processes of 16 is made. Markov
  %transition matrix is constructed with each state corresponding to a
  %(quelength, transmission number) tuple, which is made into a continuous
  %number from 0 corresponding to (0,0) and 4*(16+1) corresponding to
  %(16,4). Steady state probability is the eigen vector corresponding to
  %unit eigenvalue of the transpose of the transition probability matrix.

if(length(pVector) ==1)
    p(1:4)=pVector;
elseif(length(pVector) ==4)
    p=pVector;
else
    error("incompatible");
end

L = (Qmax+1)*M;
p_prime = 1 - p;  %complement of p
f_prime = 1 - f;  %complement of f
CSO_prime = 1 - CSO; %complement of CSO

P = zeros(L,L);
for i = 2:M
    P(i,i) = 1; %(0,m), m>= 2 unreachable states.
end

%Transition (0,1) -> (0,1)
q = 0; m = 1; 
s1 = q * M + m;
P(s1, s1) = 1 - CSO*f -CSO_prime * f * p(1);
%Transition (0,1) -> (1,1)
q = 1; m = 1; s2 = q * M + m;
P(s1, s2) = CSO*f;

%Transition (0,1) -> (1,2)
q = 1; m = 2; s2 = q * M + m;
P(s1, s2) = CSO_prime * f * p(1);

%Iterate through all values of q and m
for q = 1:Qmax
    for m = 1:M
        s1 = q * M + m;  %Current state (q,m)

        if q < Qmax && m < M
            %Next state (q,m+1)
            s2 = q * M + (m + 1);  
            P(s1, s2) =  CSO_prime*f_prime * p(m);  

            %Next state (q+1,m+1)
            s2 = (q + 1) * M + (m + 1);
            P(s1, s2) = CSO_prime*f * p(m);  
        end
        if q < Qmax
            %Next state (q+1, m) additional status for CSO
            s2 = (q + 1) * M + m;
            P(s1, s2) = CSO*f;  %some status can be (q,1)

            %Next state (q, m) additional status for CSO
            s2 = q * M + m; s2 = s1;
            P(s1, s2) = P(s1, s2) + CSO*f_prime;  %some status can be (q,1)  
        end
        if  m < M
            %Next state (q,1)
            s2 = q * M + 1;   
            P(s1, s2) = P(s1, s2)+ CSO_prime*f * p_prime(m);  

            %Next state (q-1,1)
            s2 = (q - 1) * M + 1;   
            P(s1, s2) =  CSO_prime*f_prime * p_prime(m);  

        end

        %current state (q, M)
        if m == M
            %Next state (q,1)
            s2 = q * M + 1;
            P(s1, s2) = CSO_prime*f;  

            %Next state (q-1,1)
            s2 = (q - 1) * M + 1; 
            P(s1, s2) = CSO_prime*f_prime;  
        end
        %current state (Q_max, m)
        if q == Qmax && m < M
            %Next state (Qmax,m+1)
            s2 = q * M + (m + 1);
            P(s1, s2) = CSO_prime * p(m); 
        end
        if q == Qmax
            %Next state (Q_max, m)
            P(s1, s1) = 1 - sum(P(s1, :));
        end
    end
end
[ V , D ] = eig( P' ); %Find eigenvector/eigenvalue
for i =1:size(D,1)
    %Double confirming and choosing  the eigenvector of unit eigenvalue
    if (D(i,i)>1-sqrt(eps) && D(i,i)<1+sqrt(eps)) 
        vec = V(:,i); 
        break;
    end
end

%Normalise to 1 to represent probability
vec = vec/sum(vec);

%Marginalise over m for each q to get probability of queuelengths
QLenPmf=sum(reshape(vec,M,Qmax+1));

% prev_val = 0;
% for q = 1:length(QLenPmf)
%     QLenPmf(q) = (prev_val+QLenPmf(q)) * (1-CSO);
%     prev_val = QLenPmf(q)*CSO;
% end
% % QLenPmf
end










% 
% 
% 
% 
% 
% 
% %%%Savings%%%
% L = (Qmax+1)*M;
% p_prime = 1 - p;  complement of p
% f_prime = 1 - f;  complement of f
% CSO_prime = 1 - CSO; complement of CSO
% 
% P = zeros(L,L);
% for i = 2:M
%     P(i,i) = 1; (0,m), m>=2l unreachable states.
% end
% 
% Transition (0,1) -> (0,1)
% q = 0; m = 1; 
% s1 = q * M + m;
% P(s1, s1) = 1 - CSO*f -CSO_prime * f * p(1);
% Transition (0,1) -> (1,1)
% q = 1; m = 1; s2 = q * M + m;
% P(s1, s2) = CSO*f;
% 
% Transition (0,1) -> (1,2)
% q = 1; m = 2; s2 = q * M + m;
% P(s1, s2) = CSO_prime * f * p(1);
% 
% Iterate through all values of q and m
% for q = 1:Qmax
%     for m = 1:M
%         s1 = q * M + m;  Current state (q,m)
% 
%         if q < Qmax && m < M
%             Next state (q,m+1)
%             s2 = q * M + (m + 1);  
%             P(s1, s2) = P(s1, s2)+ CSO_prime*f_prime * p(m);  
% 
%             Next state (q+1,m+1)
%             s2 = (q + 1) * M + (m + 1);
%             P(s1, s2) = P(s1, s2)+ CSO_prime*f * p(m);  
%         end
%         if q < Qmax
%             %Next state (q+1, m) additional status for CSO
%             s2 = (q + 1) * M + m;
%             P(s1, s2) = P(s1, s2)+CSO*f;  some status can be (q,1)
% 
%             %Next state (q, m) additional status for CSO
%             s2 = q * M + m; s2 = s1
%             P(s1, s2) = P(s1, s2) + CSO*f_prime;  some status can be (q,1)  
%         end
%         if  m < M
%             Next state (q,1)
%             s2 = q * M + 1;   
%             P(s1, s2) = P(s1, s2)+ CSO_prime*f * p_prime(m);  
% 
%             Next state (q-1,1)
%             s2 = (q - 1) * M + 1;   
%             P(s1, s2) = P(s1, s2)+ CSO_prime*f_prime * p_prime(m);  
% 
%         end
% 
%         current state (q, M)
%         if m == M
%             Next state (q,1)
%             s2 = q * M + 1;
%             P(s1, s2) = P(s1, s2)+CSO_prime*f;  
% 
%             Next state (q-1,1)
%             s2 = (q - 1) * M + 1; 
%             P(s1, s2) = P(s1, s2)+CSO_prime*f_prime;  
%         end
%         current state (Q_max, m)
%         if q == Qmax && m < M
%             Next state (Qmax,m+1)
%             s2 = q * M + (m + 1);
%             P(s1, s2) = P(s1, s2)+CSO_prime * p(m); 
%         end
%         if q == Qmax
%             Next state (Q_max, m)
%             P(s1, s1) = 1 - sum(P(s1, :));
%         end
%     end
% end
% [ V , D ] = eig( P' ); Find eigenvector/eigenvalue
% for i =1:size(D,1)
%     Double confirming and choosing  the eigenvector of unit eigenvalue
%     if (D(i,i)>1-sqrt(eps) && D(i,i)<1+sqrt(eps)) 
%         vec = V(:,i); 
%         break;
%     end
% end