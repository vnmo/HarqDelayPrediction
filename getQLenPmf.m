function [QLenPmf] = getQLenPmf(pVector,f,M,Qmax)

if(length(pVector) ==1)
    p(1:4)=pVector;
elseif(length(pVector) ==4)
    p=pVector;
else
    error("incompatible");
end




L = (Qmax+1)*M;
p_prime = 1 - p;  % complement of p
f_prime = 1 - f;  % complement of f

P = zeros(L,L);
for i = 2:M
    P(i,i) = 1; % (0,m), m>=2l unreachable states.
end

q = 0; m = 1; s1 = q * M + m;
P(s1, s1) = 1 - f * p(1);

q = 1; m = 2; s2 = q * M + m;
P(s1, s2) = f * p(1);

for q = 1:Qmax
    for m = 1:M
        s1 = q * M + m;  % Current state (q,m)

        if q < Qmax && m < M
            s2 = q * M + (m + 1);  
            P(s1, s2) = f_prime * p(m);  

            s2 = (q + 1) * M + (m + 1);
            P(s1, s2) = f * p(m);  
        end

        if m < M
            s2 = q * M + 1;   
            P(s1, s2) = f * p_prime(m);  

            s2 = (q - 1) * M + 1;   
            P(s1, s2) = f_prime * p_prime(m);  
        end

        if m == M
            s2 = q * M + 1;
            P(s1, s2) = f;  

            s2 = (q - 1) * M + 1; 
            P(s1, s2) = f_prime;  
        end

        if q == Qmax && m < M
            s2 = q * M + (m + 1);
            P(s1, s2) = p(m); 
        end
    end
end
P

[ V , D ] = eig( P' ); % Find eigenvector/eigenvalue

for i =1:size(D,1)
    if (D(i,i)>1-sqrt(eps) && D(i,i)<1+sqrt(eps)) 
        vec = V(:,i); 
        break;
    end
end

vec = vec/sum(vec);

QLenPmf=sum(reshape(vec,M,Qmax+1));
QLenPmf
end