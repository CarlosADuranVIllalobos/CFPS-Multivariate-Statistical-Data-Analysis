
%%Pred_accuracy.m
% Script to compare the prediction accuracy of different alg.
% With very few modifications (removing the reading process of a dataset) 
% it can calculate the accuracy of different data-sets
% The predictions change substantially to the published due to yield normalization
% for confidenciality protection
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE

%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub

clear all;
%% Initialize variables
X=[]; %Matrix of PCA variables
Y=[]; %Response variables

%% Read process measurements D1
% yield=xlsread('D1',1,'M2:M25'); %Read files 
% yield(16)=[];yield(21)=[];  %Remove not measured observations
% rtime= xlsread('D1',1,'L2:L25'); 
% rtime(16)=[];rtime(21)=[];
% DO = xlsread('D1',2,'B2:Y481')';
% DO(16,:)=[];DO(21,:)=[];
% pH = xlsread('D1',3,'B2:Y481')'; 
% pH(16,:)=[];pH(21,:)=[];
% T = xlsread('D1',4,'B2:Y481')'; 
% T(16,:)=[];T(21,:)=[];
% % Averaging every hour
% time=8; 
% sampling=480;
% step=floor(sampling/time);
% DO2=[];
% for i=1:step:sampling-step+1
%     DO2=[DO2 mean(DO(:,i:i+step-1),2) ];  
% end
% pH2=[];
% for i=1:step:sampling-step+1
%     pH2=[pH2 mean(pH(:,i:i+step-1),2) ];   
% end
% T2=[];
% for i=1:step:sampling-step+1
%     T2=[T2 mean(T(:,i:i+step-1),2) ];   
% end
% % Find index for time
% ind4=find(rtime==4); 
% ind6=find(rtime==6);
% ind8=find(rtime==8);
% % Assign values to PLS matrices
% X=[pH2(ind4,1) DO2(ind4,1) T2(ind4,1) pH2(ind4,3) DO2(ind4,4) T2(ind4,4) rtime(ind4) ];
% X=[X; pH2(ind6,1) DO2(ind6,1) T2(ind6,1) pH2(ind6,5) DO2(ind6,6) T2(ind6,6) rtime(ind6) ];  
% X=[X; pH2(ind8,1) DO2(ind8,1) T2(ind8,1) pH2(ind8,7) DO2(ind8,8) T2(ind8,8) rtime(ind8) ];
% Y=[yield(ind4);yield(ind6);yield(ind8)];

%% Read process measurements D2
% yield=xlsread('D2',1,'M2:M25'); %Read files 
% rtime= xlsread('D2',1,'L2:L25'); 
% DO = xlsread('D2',2,'B2:Y481')';
% pH = xlsread('D2',3,'B2:Y481')'; 
% T = xlsread('D2',4,'B2:Y481')'; 
% % Averaging every hour
% time=6; 
% sampling=360;
% step=floor(sampling/time);
% DO2=[];
% for i=1:step:sampling-step+1
%     DO2=[DO2 mean(DO(:,i:i+step-1),2) ];   
% end
% pH2=[];
% for i=1:step:sampling-step+1
%     pH2=[pH2 mean(pH(:,i:i+step-1),2) ];    
% end
% T2=[];
% for i=1:step:sampling-step+1
%     T2=[T2 mean(T(:,i:i+step-1),2) ];    
% end
% % Find index for length of reaction
% ind2=find(rtime==2);
% ind4=find(rtime==4);
% ind6=find(rtime==6);
% % Assign values to PCA matrix
% X=[X; pH2(ind2,1) DO2(ind2,1) T2(ind2,1) pH2(ind2,2) DO2(ind2,2) T2(ind2,2) rtime(ind2) ]; 
% X=[X; pH2(ind4,1) DO2(ind4,1) T2(ind4,1) pH2(ind4,4) DO2(ind4,4) T2(ind4,4) rtime(ind4) ];
% X=[X; pH2(ind6,1) DO2(ind6,1) T2(ind6,1) pH2(ind6,6) DO2(ind6,6) T2(ind6,6) rtime(ind6) ];
% Y=[Y; yield(ind2);yield(ind4);yield(ind6)];

%% Read process measurements D3
yield=xlsread('D3',1,'M2:M25'); %Read files 
rtime= xlsread('D3',1,'L2:L25'); 
DO = xlsread('D3',2,'B2:Y490')';
pH = xlsread('D3',3,'B2:Y490')'; 
T = xlsread('D3',4,'B2:Y490')'; 
% Averaging every hour
time=8; 
sampling=480;
step=floor(sampling/time);
DO2=[];
pH2=[];
T2=[];
for i=1:step:sampling-step+1
    DO2=[DO2 mean(DO(:,i:i+step-1),2) ];   
    pH2=[pH2 mean(pH(:,i:i+step-1),2) ]; 
    T2=[T2 mean(T(:,i:i+step-1),2) ]; 
end
% Assign values to PLS matrices
X=[X; pH2(:,1) DO2(:,1) T2(:,1) pH2(:,8) DO2(:,8) T2(:,8) rtime(:) ];
Y=[Y; yield(:)];

%% Read process measurements D4
yield=xlsread('D4',1,'M2:M25'); %Read files 
yield(4)=[];                   %Remove not measured opbservations
rtime= xlsread('D4',1,'L2:L25'); 
rtime(4)=[];  
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
DO2=[];
pH2=[];
T2=[];
for i=1:step:sampling-step+1
    DO2=[DO2 mean(DO(:,i:i+step-1),2) ];   
    pH2=[pH2 mean(pH(:,i:i+step-1),2) ]; 
    T2=[T2 mean(T(:,i:i+step-1),2) ]; 
end% Assign values to LS matrices
X=[X; pH2(:,1) DO2(:,1) T2(:,1) pH2(:,12) DO2(:,12) T2(:,12) rtime(:)];
Y=[Y; yield(:)];

 %% Read process measurements D5 (Uncomment to add)
% yield = xlsread('D5',1,'G2:G10'); %Read files 1st part 
% yield = [yield; xlsread('D5',1,'G14:G24')]; %Read files 2nd part 
% rtime= xlsread('D5',1,'E2:E10'); 
% rtime = [rtime; xlsread('D5',1,'E14:E24')];
% DO = xlsread('D5',2,'B2:J1072')'; 
% DO = [DO; xlsread('D5',2,'N2:X1072')'];
% pH = xlsread('D5',3,'B2:J1072')'; 
% pH = [pH; xlsread('D5',3,'N2:X1072')'];
% T = xlsread('D5',4,'B2:J1072')'; 
% T = [T; xlsread('D5',4,'N2:X1072')'];
% % Averaging every hour
% time=17; 
% sampling=1071;
% step=floor(sampling/time);
% DO2=[];
% for i=1:step:sampling-step+1
%     DO2=[DO2 mean(DO(:,i:i+step-1),2) ];   
% end
% pH2=[];
% for i=1:step:sampling-step+1
%     pH2=[pH2 mean(pH(:,i:i+step-1),2) ];    
% end
% T2=[];
% for i=1:step:sampling-step+1
%     T2=[T2 mean(T(:,i:i+step-1),2) ];    
% end
% % Find index for length of reaction
% ind10s=find(rtime==10);
% ind12s=find(rtime==12);
% ind14s=find(rtime==14);
% % Assign values to PLS matrices
% X=[X; pH2(ind10s,1) DO2(ind10s,1) T2(ind10s,1) pH2(ind10s,10) DO2(ind10s,10) T2(ind10s,10) rtime(ind10s) ];
% X=[X; pH2(ind12s,1) DO2(ind12s,1) T2(ind12s,1) pH2(ind12s,12) DO2(ind12s,12) T2(ind12s,12) rtime(ind12s) ]; 
% X=[X; pH2(ind14s,1) DO2(ind14s,1) T2(ind14s,1) pH2(ind14s,14) DO2(ind14s,14) T2(ind14s,14) rtime(ind14s) ];
% Y=[Y; yield(ind10s);yield(ind12s);yield(ind14s)];

%% calculate prediction errors
[X2,xmean,xstd]=zscore(X);
[Y2,ymean,ystd]=zscore(Y);
cv=plscv(X2,Y2,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
batches=length(Y);
%Calculate predictions removing the observation from the identification set
yest=[];
yest2=[];
yest3=[];
for i=1:1:batches
    yest=[yest; predm(X,Y,i,batches,lv) ]; %PLS
    yest2=[yest2; predpoly(X,Y,i,lv,2) ]; %QPLS
    yest3=[yest3; predlinear(X,Y,i,batches) ];%OLS
    
end
%Calculate SMAPE
msemlr=100*sum(abs(Y-yest3)./(0.5*(abs(Y)+abs(yest3))))/batches; 
msepls=100*sum(abs(Y-yest)./(0.5*(abs(Y)+abs(yest))))/batches;
mseqpls=100*sum(abs(Y-yest2)./(0.5*(abs(Y)+abs(yest2))))/batches;

mse=[msemlr; msepls; mseqpls];
figure(1)
bar(mse')
hold on;
xticklabels({'OLS','PLS','QPLS'})
ylabel('RMSE of prediction %');
xlabel('Length of reaction');
set(findall(gcf,'-property','FontSize'),'FontSize',34,'FontWeight','bold')

figure(2)
plot(Y,'Marker','+','Markersize',10,'linewidth',2,'LineStyle',':')
hold on
plot(yest,'Marker','o','Markersize',10,'linewidth',2,'LineStyle',':')
plot(yest2,'Marker','square','Markersize',10,'linewidth',2,'LineStyle',':')
plot(yest3,'Marker','diamond','Markersize',10,'linewidth',2,'LineStyle',':')
ylabel('Yield(normalized)')
xlabel('Runs')
legend('Measured','PLS','QPLS','OLS')
set(findall(gcf,'-property','FontSize'),'FontSize',30)
  

