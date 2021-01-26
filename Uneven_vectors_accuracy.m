
%% Uneven_vectors_accuracy.m
% Script to compare the prediction accuracy of D3 & D4 using uneven vector techniques.
% With a few modifications (removing the reading process of a dataset) 
% it can calculate the accuracy of different data-sets.
% The predictions change substantially to the published due to yield normalization
% for confidenciality protection. 

%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE

clear all;
%% Initialize variables
X=[]; %Matrix of PCA variables
Y=[]; %Response variables


%% Read process measurements D3
yield3=xlsread('D3',1,'M2:M25'); %Read files 
rtime3= xlsread('D3',1,'L2:L25'); 
DO = xlsread('D3',2,'B2:Y490')';
pH = xlsread('D3',3,'B2:Y490')'; 
T = xlsread('D3',4,'B2:Y490')'; 
% Averaging every hour
time=8; 
sampling=480;
step=floor(sampling/time);
DO23=[];
pH23=[];
T23=[];
for i=1:step:sampling-step+1
    DO23=[DO23 mean(DO(:,i:i+step-1),2) ];
    pH23=[pH23 mean(pH(:,i:i+step-1),2) ]; 
    T23=[T23 mean(T(:,i:i+step-1),2) ];  
end

%% Read process measurements D4
yield4=xlsread('D4',1,'M2:M25'); %Read files 
yield4(4)=[];                   %Remove not measured opbservations
rtime4= xlsread('D4',1,'L2:L25'); 
rtime4(4)=[];  
DO = xlsread('D4',2,'B2:Y721')';
DO(4,:)=[];
pH = xlsread('D4',3,'B2:Y721')';
pH(4,:)=[];
T = xlsread('D4',4,'B2:Y721')'; 
T(4,:)=[];
% Averaging every hour
time=12; 
sampling=720;
step=floor(sampling/time)-1;
DO24=[];
pH24=[];
T24=[];
for i=1:step:sampling-step+1
    DO24=[DO24 mean(DO(:,i:i+step-1),2) ];
    pH24=[pH24 mean(pH(:,i:i+step-1),2) ]; 
    T24=[T24 mean(T(:,i:i+step-1),2) ];  
end
yield=[yield3;yield4];
rtime=[rtime3;rtime4];
%% Preditions using fist 8h (Approach 1)
time=8;
X=[];
for i=1:1:length(yield3)
Xn=[];
    for n=1:1:time
        Xn = [Xn  pH23(i,n) DO23(i,n) T23(i,n)];
    end
     X=[X;Xn];
end
for i=1:1:length(yield4)
Xn=[];
    for n=1:1:time
        Xn = [Xn  pH24(i,n) DO24(i,n) T24(i,n)];
    end
    X=[X;Xn];
end
X=[X, rtime];
batches=length(yield);
Y=[];
for i=1:1:batches
     Yn=[];
     Yn = [Yn yield(i)];
     Y=[Y;Yn];
end
[X2,xmean,xstd]=zscore(X);
[Y2,ymean,ystd]=zscore(Y);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
%calculate predictions
yest=[];
yest2=[];
yest3=[];
for i=1:1:batches
    yest=[yest; predm(X,Y,i,batches,lv) ]; %PLS
    yest2=[yest2; predpoly(X,Y,i,lv,2) ]; %QPLS
    yest3=[yest3; predlinear(X,Y,i,batches) ];%OLS
end
%calculate SMAPE
msemlrf=100*sum(abs(Y-yest3)./(0.5*(abs(Y)+abs(yest3))))/batches; 
mseplsf=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
mseqplsf=100*sum(abs(Y-yest2)./(0.5*(abs(Y)+abs(yest2))))/batches;

%% Preditions using last 8h (Approach 2)
time=8;
X=[];
for i=1:1:length(yield3)
Xn=[];
    for n=1:1:time
        Xn = [Xn  pH23(i,n) DO23(i,n) T23(i,n)];
    end
     X=[X;Xn];
end
time=12;
for i=1:1:length(yield4)
Xn=[];
    for n=5:1:time
        Xn = [Xn  pH24(i,n) DO24(i,n) T24(i,n)];
    end
    X=[X;Xn];
end
X=[X, rtime];
Y=[];
for i=1:1:batches
     Yn=[];
     Yn = [Yn yield(i)];
     Y=[Y;Yn];
end
[X2,xmean,xstd]=zscore(X);
[Y2,ymean,ystd]=zscore(Y);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
%calculate predictions
yest=[];
yest2=[];
yest3=[];
for i=1:1:batches
    yest=[yest; predm(X,Y,i,batches,lv) ]; %PLS
    yest2=[yest2; predpoly(X,Y,i,lv,2) ]; %QPLS
    yest3=[yest3; predlinear(X,Y,i,batches) ];%OLS
end
%calculate SMAPE
msemlrl=100*sum(abs(Y-yest3)./(0.5*(abs(Y)+abs(yest3))))/batches; 
mseplsl=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
mseqplsl=100*sum(abs(Y-yest2)./(0.5*(abs(Y)+abs(yest2))))/batches;

%% Preditions using missing data techniques (Approach 3)
time=12;
X=[];
for i=1:1:23
Xn=[];
    for n=1:1:time
        Xn = [Xn  pH24(i,n) DO24(i,n) T24(i,n)];
    end
     X=[X;Xn];
end
[X2,xmean,xstd]=zscore(X);
%calculate future values for D3
[P,T,latent,tsquared,explained,mu] = pca(X2);
pc=3; %Number of PC from explained 
Xnew=[];
time=8;
for i=1:1:24
Xn=[];
    for n=1:1:time
        Xn = [Xn  pH23(i,n) DO23(i,n) T23(i,n)];
    end
     Xnew=[Xnew;Xn];
end
Xnew2=(Xnew-xmean(1:24))./xstd(1:24);
Tnew=(P(1:24,1:3)'*P(1:24,1:3))\P(1:24,1:3)'*Xnew2';
Xnewfuture=Tnew'*P(25:end,1:3)';
Xnewf=Xnewfuture.*xstd(25:end)+xmean(25:end);
Xn2=[Xnew Xnewf]; %Add future values to D3
X=[X;Xn2];
X=[X,rtime];
[X2,xmean,xstd]=zscore(X);
[Y2,ymean,ystd]=zscore(Y);

cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
%calculate predictions
yest=[];
yest2=[];
yest3=[];
for i=1:1:batches
    yest=[yest; predm(X,Y,i,batches,lv) ]; %PLS
    yest2=[yest2; predpoly(X,Y,i,lv,2) ]; %QPLS
    yest3=[yest3; predlinear(X,Y,i,batches) ];%OLS
end
%calculate SMAPE
msemlrm=100*sum(abs(Y-yest3)./(0.5*(abs(Y)+abs(yest3))))/batches; 
mseplsm=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
mseqplsm=100*sum(abs(Y-yest2)./(0.5*(abs(Y)+abs(yest2))))/batches;

%% Plot accurace errors
mse=[msemlrf mseplsf mseqplsf ; msemlrl mseplsl mseqplsl ; msemlrm mseplsm mseqplsm ];   
     
 
figure(5)
    bar(mse)
  hold on;
   xticklabels({'First 8h','Last 8h','Missing data'})
   xtickangle(45)
   ylabel('RMSE of prediction %');
   legend({'MLR','PLS','QPLS'})
   grid on
  set(findall(gcf,'-property','FontSize'),'FontSize',28)  

