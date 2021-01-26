function yest=predpoly(X,Y,stop,lv,order)
% Function to do leave-one-out predictions using Polynomial regression
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE
%   yest: Predicted responses
%   X: Predictor variables
%   yield: Response variables
%   stop: leave-one-out response index
%   batches: number of response observations
%   lv: Number of Latent Variables
%   order:Order of the polynomial
xk=X(stop,:);

X(stop,:)=[];
Y(stop,:)=[];


[m,n] = size(X);
mxcal    = mean(X);
mcxcal   = (X-mxcal(ones(m,1),:));

[m,n] = size(Y);
mycal    = mean(Y);
mcycal   = (Y-mycal(ones(m,1),:));


[p,q,w,t,u,b2] = polypls(mcxcal,mcycal,lv,order);
[m,n] = size(xk);
scxtest  = xk-mxcal(ones(m,1),:);
yest = polypred(scxtest,b2,p,q,w,lv);
[m,n] = size(yest);
yest  =yest+mycal(ones(m,1),:);
end


