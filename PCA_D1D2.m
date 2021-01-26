
%%PCA_D1D2.m
% Script to do PCA analysis of a series of Cell-Free Synthesis (D1 and D2)

%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% PCA Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE
clear all;
%% Initialize variables
X=[]; %Matrix of PCA variables

categories=['pH ini'     %Variable labels for the plots 
            'DO ini' 
            'T  ini'   
            'pH end'    
            'DO end' 
            'T  end' 
            'Length'
            'Yield '     
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
% Assign values to PCA matrix
X=[pH2(ind4,1) DO2(ind4,1) T2(ind4,1) pH2(ind4,3) DO2(ind4,4) T2(ind4,4) rtime(ind4) yield(ind4)];
X=[X; pH2(ind6,1) DO2(ind6,1) T2(ind6,1) pH2(ind6,5) DO2(ind6,6) T2(ind6,6) rtime(ind6) yield(ind6)];  
X=[X; pH2(ind8,1) DO2(ind8,1) T2(ind8,1) pH2(ind8,7) DO2(ind8,8) T2(ind8,8) rtime(ind8) yield(ind8)];

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
% Assign values to PCA matrix
X=[X; pH2(ind2,1) DO2(ind2,1) T2(ind2,1) pH2(ind2,2) DO2(ind2,2) T2(ind2,2) rtime(ind2) yield(ind2)]; 
X=[X; pH2(ind4,1) DO2(ind4,1) T2(ind4,1) pH2(ind4,4) DO2(ind4,4) T2(ind4,4) rtime(ind4) yield(ind4)];
X=[X; pH2(ind6,1) DO2(ind6,1) T2(ind6,1) pH2(ind6,6) DO2(ind6,6) T2(ind6,6) rtime(ind6) yield(ind6)];


%% Obtain a PCA model
[X2,xmean,xstd]=zscore(X);%Normalize
[wcoeff,score,latent,tsquared,explained,mu] = pca(X2);
% Plot the variance explained by each PC
figure(1)
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
set(findall(gcf,'-property','FontSize'),'FontSize',30)

%% Plot by experiment
figure(2)
h=biplot(wcoeff(:,1:2),'scores',score(:,1:2),'varlabels',categories,'LineWidth',1.5,'marker','s');
hID = get(h, 'tag');
hPt = h(strcmp(hID,'obsmarker'));
%D1 
NameArray = {'Marker','color','Markersize','MarkerFaceColor','DisplayName'};
ValueArray = {'o','blue',6,'blue','D1 4-8h'};
set(hPt(1:22),NameArray,ValueArray);
%D2
NameArray = {'Marker','color','Markersize','MarkerFaceColor','DisplayName'};
ValueArray = {'diamond','magenta',6,'magenta','D2 2-6h'};
set(hPt(23:46),NameArray,ValueArray);

legend([hPt(1),hPt(23)])
xlabel('PC 1')
ylabel('PC 2')

set(findall(gcf,'-property','FontSize'),'FontSize',28)
box on;

