function [U,F] = solve_YACHT(H_ini,numanchor,w,YY_label)
maxIter = 50 ; 
[n,cc] = size(H_ini); 
m = numanchor;
res_base = numanchor;
max_mu = 1e8;
A = zeros(res_base,m);         % m  * m
Z = zeros(m,n);          % m  * n
P = ones(m,cc);
rho = 1.3;  %adjust \rho to obatain optimal performance. rho =  1.1 1.2 1.4
mu = 0.5;  %adjust \mu to obatain optimal performance. mu = 1 mu = 0.1
one_m = ones(m,1);
one_n = ones(n,1);
J = ones(n,1);
YY = YY_label;  %

flag = 1;
iter = 0;
F = orth(YY);
rho_1=1;
H = H_ini';
clear  YY_label X
while flag
    
    iter = iter + 1; 

    %% optimize P
    
    Pt = P;
    AZ = A*Z; % since each view share the same A and Z.
     
    C = AZ*H';      
    
    [Up,~,Vp] = svd(C,'econ');
    P = Up*Vp';%
    %% optimize F

    Ft = F;
    [F] = updateF(H',F,YY,w);  


    %% optimize A
    At = A;
    PX = P*H; % since each view share the same A and Z.
    G = PX*Z';      
    [Ug,~,Vg] = svd(G,'econ');
    A = Ug*Vg';

    %% optimize Z
    Zt = Z;
    temp = 2*A'*P*P'*A + mu.* one_m*one_m' ;
    Z = (temp) \ (mu.*(one_m*one_n')+2*A'*P*H-one_m*J');


    %% optimize H
    Ht = H;
    H = gradient_descent(H,w,F,rho_1,P,A,Z);

    
    %% optimize J and mu
    J = J + mu*(Z'*one_m-one_n);
    mu = min(max_mu,mu*rho);

    %% optimize alpha

    diffC = abs(norm(P - Pt,'fro')/norm(Pt,'fro'));
    diffE = abs(norm(A - At,'fro')/norm(At,'fro'));
    diffF = abs(norm(F - Ft,'fro')/norm(Ft,'fro'));
    diffH = abs(norm(H - Ht,'fro')/norm(Ht,'fro'));
    diffJ = abs(norm(Z - Zt,'fro')/norm(Zt,'fro'));
    stopC = max([diffC,diffE,diffF,diffJ,diffH]);
    if iter>2
    if  (stopC<1e-1 || iter>maxIter)
        [U,~,~]=svd(Z','econ');
        break
    end
    end

end

end         
