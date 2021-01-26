function yest=predlinear(X,yield,stop,batches)
% Function to do leave-one-out predictions using OLS regression
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


Y=[];
Xs=[];
Y=[];
Xs=[];

for i=1:1:batches
    if(i==stop)
    else
     Yn=[];
     Yn = [Yn yield(i)];
     
     Y=[Y;Yn]; 
     Xs=[Xs ; X(i,:)];
    end
    
end

[X2,xmean,xstd]=zscore(Xs);
[Y2,ymean,ystd]=zscore(Y);
 xstd(xstd<=1e-12)=1e-12;
 
 mdl = fitlm(X2,Y2);
Xp=(X(stop,:)-xmean)./xstd;
yp = predict(mdl,Xp);

yest= yp.*ystd+ymean;
end

%+++ END ++++++++++++++++++++++++++++++++++++
%+++ There is a song you like to sing in your memory.
