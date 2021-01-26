function [B,Wstar,T,P,Q,W,R2X,R2Y]=plsnipals(X,Y,lv)
% Function to predict using PLS regression using the NIPALS algorithm.
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE
%   B: Coefficiets vector
%   Wstar: x-block weights to calculate new scores
%   P: x-block loadings
%   Q: y-block loadings
%   W: x-block weights
%   T: x-block scores
%   R2X: explained variance in X
%   R2Y: explained variance in Y
%   Y: Matrix of the predicted variables
%   X: matrix of predictor variables
%   lv: Number of Latent Variables




varX=sum(sum(X.^2));
varY=sum(sum(Y.^2));
for i=1:lv
    err=1;
    u=Y(:,1);
    niter=0;
    while (err>1e-8 && niter<1000)  % for convergence 
        w=X'*u/(u'*u);
        w=w/norm(w);
        t=X*w;
        q=Y'*t/(t'*t);              
        u1=Y*q/(q'*q);
        err=norm(u1-u)/norm(u);
        u=u1;
        niter=niter+1;
    end
    p=X'*t/(t'*t);
    X=X-t*p';
    Y=Y-t*q';
    
    W(:,i)=w;
    T(:,i)=t;
    P(:,i)=p;
    Q(:,i)=q;
    
end

% calculate explained variance
R2X=diag(T'*T*P'*P)/varX;
R2Y=diag(T'*T*Q'*Q)/varY;
%calculate the transformed W matrix
Wstar=W*(P'*W)^(-1); 
B=Wstar*Q';
Q=Q';

%+++ 
