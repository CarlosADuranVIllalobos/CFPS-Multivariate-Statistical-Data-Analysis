
%%Pred_accuracy_bytime.m
% Script to compare the prediction accuracy of D5 using different time slots.
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


 %% Read process measurements D5 
yield = xlsread('D5',1,'G2:G10'); %Read files 1st part 
yield = [yield; xlsread('D5',1,'G14:G24')]; %Read files 2nd part 
mon = xlsread('D5',1,'H2:H10'); %Read files 1st part 
mon = [mon; xlsread('D5',1,'H14:H24')]; %Read files 2nd part 
agg = xlsread('D5',1,'I2:I10'); %Read files 1st part 
agg = [agg; xlsread('D5',1,'I14:I24')]; %Read files 2nd part 
rtime= xlsread('D5',1,'E2:E10'); 
rtime = [rtime; xlsread('D5',1,'E14:E24')];
DO = xlsread('D5',2,'B2:J1072')'; 
DO = [DO; xlsread('D5',2,'N2:X1072')'];
pH = xlsread('D5',3,'B2:J1072')'; 
pH = [pH; xlsread('D5',3,'N2:X1072')'];
T = xlsread('D5',4,'B2:J1072')'; 
T = [T; xlsread('D5',4,'N2:X1072')'];
O2=xlsread('D5',5,'B2:J1072')'; 
O2 = [O2; xlsread('D5',5,'N2:X1072')'];
CO2=xlsread('D5',6,'B2:J1072')'; 
CO2 = [CO2; xlsread('D5',6,'N2:X1072')'];
NH3=xlsread('D5',7,'B2:J1072')'; 
NH3 = [NH3; xlsread('D5',7,'N2:X1072')'];

% Averaging every hour
time=17; 
sampling=1071;
step=floor(sampling/time);
DO2=[];
pH2=[];
T2=[];
O22=[];
CO22=[];
NH32=[];
for i=1:step:sampling-step+1
    DO2=[DO2 mean(DO(:,i:i+step-1),2) ]; 
    pH2=[pH2 mean(pH(:,i:i+step-1),2) ];
    T2=[T2 mean(T(:,i:i+step-1),2) ];
    O22=[O22 mean(O2(:,i:i+step-1),2) ];
    CO22=[CO22 mean(CO2(:,i:i+step-1),2) ];  
    NH32=[NH32 mean(NH3(:,i:i+step-1),2) ];  
end

% Assign values to PLS matrices
batches=length(yield);
time=14; %number of hours
X=[];
Y=[];
for i=1:1:batches
Xn=[];
   for n=1:1:time
        Xn = [Xn  pH2(i,n) DO2(i,n) T2(i,n) O22(i,n) NH32(i,n)];     
        if n>=5
         Xn = [Xn CO22(i,n)]; %assign real values only 
        end
    end
X=[X;Xn];
Yn=[];  
Yn = [Yn mon(i)]; %Change to yield or agg. for different response variables
Y=[Y;Yn];
end

% calculate prediction errors
%% Using 1h
Xn=X(:,1:5);
[X2,xmean,xstd]=zscore(Xn);
[Y2,ymean,ystd]=zscore(Y);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls1=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 2h  
Xn=X(:,1:10);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls2=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 3h  
Xn=X(:,1:15);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls3=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 4h  
Xn=X(:,1:20);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls4=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 5h  
Xn=X(:,1:26);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls5=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 6h  
Xn=X(:,1:32);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls6=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 7h  
Xn=X(:,1:38);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls7=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 8h  
Xn=X(:,1:44);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls8=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 9h  
Xn=X(:,1:50);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls9=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 10h  
Xn=X(:,1:56);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls10=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 11h  
Xn=X(:,1:62);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls11=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
%% Using 12h  
Xn=X(:,1:68);
[X2,xmean,xstd]=zscore(Xn);
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
yest=[];
for i=1:1:batches
    yest=[yest; predm(Xn,Y,i,batches,lv) ];
end
yieldmsepls12=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;

%% Plot results
mse=[yieldmsepls1 yieldmsepls2 yieldmsepls3 yieldmsepls4 yieldmsepls5 yieldmsepls6 yieldmsepls7 yieldmsepls8 yieldmsepls9 yieldmsepls10 yieldmsepls11 yieldmsepls12];
figure(5)
  bar(mse)
  hold on;
  title('Monomer % SMAPE using MPLS')
  % xticklabels({'1','2','1&2 Missing Data','1&2 First 6h'})
   ylabel('RMSE of prediction %');
   xlabel('Time (h)');
   %legend({'MLR','PLS','QPLS','ENCV'})
   grid on
  set(findall(gcf,'-property','FontSize'),'FontSize',28)
