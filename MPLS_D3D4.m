
%%MPLS_D3D4.m
% Script to run a Multiway-PLS regression and confidence limits calculation of a series of Cell-Free Synthesis (D3 & D4) 

%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE

clear all;
%% Initialize variables

Xd3=[]; %Predictor variables D3
Yd3=[]; %Response variables D3
Xd4=[]; %Predictor variables D4
Yd4=[]; %Response variables D4


%% Read process measurements and obtain regression vector for D3 
yield=xlsread('D3',1,'M2:M25'); %Read files 
rtimed3= xlsread('D3',1,'L2:L25'); 
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

% Assign values the matrix for MPLS
time=8;
for i=1:1:length(yield)
Xn=[];
    for n=1:1:time
        Xn = [Xn  pH2(i,n) DO2(i,n) T2(i,n)];
    end
     Xd3=[Xd3;Xn];
     Yd3=[Yd3;yield(i)];
end
% Obtain a PLS model regression vector and confidence limits
[X2,xmean,xstd]=zscore(Xd3); %Normalize
[Y2,ymean,ystd]=zscore(Yd3); %Normalize
cv=plscv(X2,Y2,10,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
[B,limit]=clim(X2,Y2,lv); %Obtain the regression vector and confidence limit
Bph=[];
sph=[];
for i=1:3:length(B)
    Bph=[Bph B(i)]; 
    sph=[sph limit(i)]; 
end
Bdo=[];
sdo=[];
for i=2:3:length(B)
    Bdo=[Bdo B(i)]; 
    sdo=[sdo limit(i)];
end  
Bt=[];
st=[];
for i=3:3:length(B)
     Bt=[Bt B(i)]; 
     st=[st limit(i)]; 
end
%Plot regression vector
figure()
subplot(2,2,1)
 bar(Bph)
 hold on
errorbar(Bph,sph,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
 title('(a)')
grid on
box on
axis([0 9 -0.2 0.2])
  xlabel('Time h');
 ylabel('D3 Yield \beta_p_H');
subplot(2,2,2)
 bar(Bdo)
 hold on
errorbar(Bdo,sdo,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
 title('(b)')
grid on
box on
axis([0 9 -0.2 0.2])
xtickangle(60)
xlabel('Time h');
 ylabel('D3 Yield \beta_D_O');
subplot(2,2,3)
 bar(Bt)
  hold on
errorbar(Bt,st,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
 title('(c)')
grid on
box on
axis([0 9 -0.2 0.2])
xtickangle(60)
xlabel('Time h');
ylabel('D3 Yield \beta_T');
set(findall(gcf,'-property','FontSize'),'FontSize',24)


%% Read process measurements obtain regression vector for D4
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
step=floor(sampling/time);
DO2=[];
pH2=[];
T2=[];
for i=1:step:sampling-step+1
    DO2=[DO2 mean(DO(:,i:i+step-1),2) ];   
    pH2=[pH2 mean(pH(:,i:i+step-1),2) ]; 
    T2=[T2 mean(T(:,i:i+step-1),2) ]; 
end
% Assign values the matrix for MPLS
time=12;
for i=1:1:length(yield)
    Xn=[];
    for n=1:1:time
        Xn = [Xn  pH2(i,n) DO2(i,n) T2(i,n)];
    end
    
    Xd4=[Xd4;Xn];
    Yd4=[Yd4;yield(i)];
end
% Obtain a PLS model regression vector and confidence limits
[X2,xmean,xstd]=zscore(Xd4); %Normalize
[Y2,ymean,ystd]=zscore(Yd4); %Normalize
cv=plscv(X2,Y2,10,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV
[B,limit]=clim(X2,Y2,lv); %Obtain the regression vector and confidence limit
Bph=[];
sph=[];
for i=1:3:length(B)
    Bph=[Bph B(i)]; 
    sph=[sph limit(i)]; 
end
Bdo=[];
sdo=[];
for i=2:3:length(B)
    Bdo=[Bdo B(i)]; 
    sdo=[sdo limit(i)];
end  
Bt=[];
st=[];
for i=3:3:length(B)
     Bt=[Bt B(i)]; 
     st=[st limit(i)]; 
end
%Plot regression vector
figure()
subplot(2,2,1)
bar(Bph)
hold on
errorbar(Bph,sph,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
title('(a)')
grid on
box on
axis([0 13 -0.25 0.25])
xlabel('Time h');
ylabel('D3&D4 Yield \beta_p_H');

subplot(2,2,2)
bar(Bdo)
hold on
errorbar(Bdo,sdo,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
title('(b)')
grid on
box on
axis([0 13 -0.25 0.25])
xtickangle(60)
xlabel('Time h');
ylabel('D3&D4 Yield \beta_D_O');

subplot(2,2,3)
bar(Bt)
hold on
errorbar(Bt,st,'LineStyle','none','color','k','CapSize',10,'LineWidth',2)
title('(c)')
grid on
box on
axis([0 13 -0.25 0.25])
xtickangle(60)
xlabel('Time h');
ylabel('D3&D4 Yield \beta_T');
set(findall(gcf,'-property','FontSize'),'FontSize',24)


