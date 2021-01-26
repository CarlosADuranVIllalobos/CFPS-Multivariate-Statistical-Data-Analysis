
%%PLS_D1D2D5.m
% Script to run a PLS regression and confidence limits calculation of a series of Cell-Free Synthesis (D1 & D2, D5) 

%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE

clear all;
%% Initialize variables
X=[]; %Predictor variables
Y=[]; %Response variables

categoriesd1d2=['pH ini'    
            'DO ini' 
            'T  ini'   
            'pH sp '    
            'DO end' 
            'T  sp ' 
            'Length'
       ];
   
   
categoriesd5=[' pH ini   '    
            ' DO ini   ' 
            ' T  ini   '   
            ' pH sp    '    
            ' DO end   ' 
            ' T  sp    ' 
            ' Length   '
            ' O_2 v end'
            'CO_2 v end'
            'NH_3 v end'
  
       ];

%% Read process measurements D1
yield=xlsread('D1',1,'M2:M25'); %Read files 
yield(16)=[];yield(21)=[];  %Remove not measured observations
rtime= xlsread('D1',1,'L2:L25'); 
rtime(16)=[];rtime(21)=[];
DO = xlsread('D1',2,'B2:Y481')';
DO(16,:)=[];DO(21,:)=[];
pH = xlsread('D1',3,'B2:Y481')'; 
pH(16,:)=[];pH(21,:)=[];
T = xlsread('D1',4,'B2:Y481')'; 
T(16,:)=[];T(21,:)=[];
% Averaging every hour
time=8; 
sampling=480;
step=floor(sampling/time);
DO2=[];
for i=1:step:sampling-step+1
    DO2=[DO2 mean(DO(:,i:i+step-1),2) ];  
end
pH2=[];
for i=1:step:sampling-step+1
    pH2=[pH2 mean(pH(:,i:i+step-1),2) ];   
end
T2=[];
for i=1:step:sampling-step+1
    T2=[T2 mean(T(:,i:i+step-1),2) ];   
end
% Find index for time
ind4=find(rtime==4); 
ind6=find(rtime==6);
ind8=find(rtime==8);
% Assign values to the PLS matrices
X=[pH2(ind4,1) DO2(ind4,1) T2(ind4,1) pH2(ind4,3) DO2(ind4,4) T2(ind4,4) rtime(ind4)];
X=[X; pH2(ind6,1) DO2(ind6,1) T2(ind6,1) pH2(ind6,5) DO2(ind6,6) T2(ind6,6) rtime(ind6)];  
X=[X; pH2(ind8,1) DO2(ind8,1) T2(ind8,1) pH2(ind8,7) DO2(ind8,8) T2(ind8,8) rtime(ind8)];
Y=[yield(ind4);yield(ind6);yield(ind8)];

%% Read process measurements D2
yield=xlsread('D2',1,'M2:M25'); %Read files 
rtime= xlsread('D2',1,'L2:L25'); 
DO = xlsread('D2',2,'B2:Y481')';
pH = xlsread('D2',3,'B2:Y481')'; 
T = xlsread('D2',4,'B2:Y481')'; 
% Averaging every hour
time=6; 
sampling=360;
step=floor(sampling/time);
DO2=[];
for i=1:step:sampling-step+1
    DO2=[DO2 mean(DO(:,i:i+step-1),2) ];   
end
pH2=[];
for i=1:step:sampling-step+1
    pH2=[pH2 mean(pH(:,i:i+step-1),2) ];    
end
T2=[];
for i=1:step:sampling-step+1
    T2=[T2 mean(T(:,i:i+step-1),2) ];    
end
% Find index for length of reaction
ind2=find(rtime==2);
ind4=find(rtime==4);
ind6=find(rtime==6);
% Assign values to the PLS matrices
X=[X; pH2(ind2,1) DO2(ind2,1) T2(ind2,1) pH2(ind2,2) DO2(ind2,2) T2(ind2,2) rtime(ind2)]; 
X=[X; pH2(ind4,1) DO2(ind4,1) T2(ind4,1) pH2(ind4,4) DO2(ind4,4) T2(ind4,4) rtime(ind4)];
X=[X; pH2(ind6,1) DO2(ind6,1) T2(ind6,1) pH2(ind6,6) DO2(ind6,6) T2(ind6,6) rtime(ind6)];
X1= X; %Save values for D1 and D2 
Y=[Y; yield(ind2);yield(ind4);yield(ind6)]; 
Y1= Y; %Save values for D1 and D2

%% Read process measurements D5
X=[]; %Predictor variables
Y=[]; %Response variables
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
for i=1:step:sampling-step+1
    DO2=[DO2 mean(DO(:,i:i+step-1),2) ];   
end
pH2=[];
for i=1:step:sampling-step+1
    pH2=[pH2 mean(pH(:,i:i+step-1),2) ];    
end
T2=[];
for i=1:step:sampling-step+1
    T2=[T2 mean(T(:,i:i+step-1),2) ];    
end
O22=[];
for i=1:step:sampling
    O22=[O22 mean(O2(:,i:i+step-1),2) ];    
end
O2=O22;
CO22=[];
for i=1:step:sampling
    CO22=[CO22 mean(CO2(:,i:i+step-1),2) ];    
end
C02=CO22;
NH32=[];
for i=1:step:sampling
    NH32=[NH32 mean(NH3(:,i:i+step-1),2) ];  
end
NH3=NH32;

% Find index for length of reaction
ind10s=find(rtime==10);
ind12s=find(rtime==12);
ind14s=find(rtime==14);
% Assign values to PCA matrix
X=[X; pH2(ind10s,1) DO2(ind10s,1) T2(ind10s,1) pH2(ind10s,10) DO2(ind10s,10) T2(ind10s,10) rtime(ind10s) O2(ind10s,10) CO2(ind10s,10) NH3(ind10s,10) ];
X=[X; pH2(ind12s,1) DO2(ind12s,1) T2(ind12s,1) pH2(ind12s,12) DO2(ind12s,12) T2(ind12s,12) rtime(ind12s) O2(ind12s,12) CO2(ind12s,12) NH3(ind12s,12)]; 
X=[X; pH2(ind14s,1) DO2(ind14s,1) T2(ind14s,1) pH2(ind14s,14) DO2(ind14s,14) T2(ind14s,14) rtime(ind14s) O2(ind14s,14) CO2(ind14s,14) NH3(ind14s,14)];
X5= X; %Save values for D5 
Y=[Y; yield(ind10s);yield(ind12s);yield(ind14s)]; %change here for mon or agg instead of yield.
Y5= Y; %Save values for D5

%% Obtain a PLS model regression vector and confidence limits for D1 and D2
[X2,xmean,xstd]=zscore(X1); %Normalize
[Y2,ymean,ystd]=zscore(Y1); %Normalize
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
[B,limit]=clim(X2,Y2,lv); %Obtain the regression vector and confidence limit

figure(1)
subplot(1,2,1)
bar(B)
hold on
errorbar(B,limit,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
title('(a)') 
xticks([1:7]);
xticklabels(categoriesd1d2)
xtickangle(60)
axis([0 8 -0.5 0.5])
box on
grid on
yticks([-.5 -.4 -.3 -.2 -.1 0 .1 .2 .3 .4 .5])
yticklabels({'-.5','-.4','-.3','-.2','-.1','0','.1','.2','.3','.4','.5'})
ylabel('D1&D2  PLS \beta_Y_i_e_l_d');
set(findall(gcf,'-property','FontSize'),'FontSize',24)


%% Obtain a PLS model regression vector and confidence limits for D5

[X2,xmean,xstd]=zscore(X5); %Normalize
[Y2,ymean,ystd]=zscore(Y5); %Normalize
cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
[B,limit]=clim(X2,Y2,lv); %Obtain the regression vector and confidence limit

subplot(1,2,2)
bar(B)
hold on
errorbar(B,limit,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
title('(b)')
xticks([1:10]);
set(gca,'TickLabelInterpreter','tex');
xticklabels(categoriesd5);
xtickangle(60);
axis([0 11 -0.7 0.7]);
box on;
grid on;
ylabel('D5 PLS \beta_Y_i_e_l_d');
 set(findall(gcf,'-property','FontSize'),'FontSize',24);


