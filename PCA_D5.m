
%%PCA_D5.m
% Script to do PCA analysis of a series of Cell-Free Synthesis (D5)
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
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


%% Read process measurements D5
yield = xlsread('D5',1,'G2:G10'); %Read files 1st part 
yield = [yield; xlsread('D5',1,'G14:G24')]; %Read files 2nd part 
rtime= xlsread('D5',1,'E2:E10'); 
rtime = [rtime; xlsread('D5',1,'E14:E24')];
DO = xlsread('D5',2,'B2:J1072')'; 
DO = [DO; xlsread('D5',2,'N2:X1072')'];
pH = xlsread('D5',3,'B2:J1072')'; 
pH = [pH; xlsread('D5',3,'N2:X1072')'];
T = xlsread('D5',4,'B2:J1072')'; 
T = [T; xlsread('D5',4,'N2:X1072')'];
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
% Find index for length of reaction
ind10s=find(rtime==10);
ind12s=find(rtime==12);
ind14s=find(rtime==14);
% Assign values to PCA matrix
X=[X; pH2(ind10s,1) DO2(ind10s,1) T2(ind10s,1) pH2(ind10s,10) DO2(ind10s,10) T2(ind10s,10) rtime(ind10s) yield(ind10s)];
X=[X; pH2(ind12s,1) DO2(ind12s,1) T2(ind12s,1) pH2(ind12s,12) DO2(ind12s,12) T2(ind12s,12) rtime(ind12s) yield(ind12s)]; 
X=[X; pH2(ind14s,1) DO2(ind14s,1) T2(ind14s,1) pH2(ind14s,14) DO2(ind14s,14) T2(ind14s,14) rtime(ind14s) yield(ind14s)];

%% Obtain a PCA model
[X2,xmean,xstd]=zscore(X); %Normalize
[wcoeff,score,latent,tsquared,explained,mu] = pca(X2);
% Plot the variance explained by each PC
figure(1)
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')
set(findall(gcf,'-property','FontSize'),'FontSize',30)

%% Plot by experiment
figure(2)
h=biplot(wcoeff(:,1:3),'scores',score(:,1:3),'varlabels',categories,'LineWidth',1.5,'marker','s');
hID = get(h, 'tag');
hPt = h(strcmp(hID,'obsmarker'));

NameArray = {'Marker','color','Markersize','MarkerFaceColor','DisplayName'};
ValueArray = {'x',[0.5 .25 .06],6,[0.5 .25 .06],'D5 10-14H'};
set(hPt(1:end),NameArray,ValueArray);
 
legend([hPt(1)])
xlabel('PC 1')
ylabel('PC 2')

set(findall(gcf,'-property','FontSize'),'FontSize',28)
box on;

