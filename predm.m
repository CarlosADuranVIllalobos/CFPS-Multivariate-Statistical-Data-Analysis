function yest=predm(X,yield,stop,batches,lv)
% Function to do leave-one-out predictions using PLS regression
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
%   lv=Number of Latent Variables
Y=[];
Xs=[];

for i=1:1:batches
    if(i==stop)
    else
     Yn=[];
     Yn = [Yn yield(i,:)];
     
     Y=[Y;Yn]; 
     Xs=[Xs ; X(i,:)];
    end
    
end

[X2,xmean,xstd]=zscore(Xs);
[Y2,ymean,ystd]=zscore(Y);
 xstd(xstd<=1e-12)=1e-12;
 
%[P,C,T,U,beta,PCTVAR,MSE,stats] = plsregress(X2,Y2,lv);
[B,Wstar,T,P,Q,W,R2X,R2Y]=plsnipals(X2,Y2,lv);
 Xp=(X(stop,:)-xmean)./xstd;
     %yp=[1 Xp]*beta; 
     yp=Xp*B;
yest= yp.*ystd+ymean;
end

%+++ END ++++++++++++++++++++++++++++++++++++
%+++ There is a song you like to sing in your memory.
