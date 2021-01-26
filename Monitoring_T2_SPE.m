
%% Monitoring_T2_SPE.m
% Script to use PLS monitoring techniques
% D3 & D4.
% With a few modifications it can calculate the accuracy of different datasets.

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
batches=length(yield);

%% Build a dataset using last 8h (Approach 2)
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

%% Divide in 'training','good','bad' 
% These index were obtained from the non-normalized values (from the article)
% using the following rules:
% indmax=find(Y>=0.3); %For acceptable values of yield
% indmin=find(Y<=0.2); %For non-acceptable values of yield
indmax=[1;2;3;4;6;7;10;13;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47];
indmin=[5;11];
%suboptimal dataset
Xbad=(X(indmin,:)); 
Ybad=(Y(indmin));  
%training dataset
Xtr=X(indmax,:);
Ytr=Y(indmax);
%Optimal dataset
Xgood=(Xtr(21:23,:));
Ygood=(Ytr(21:23));
%remove optimal dataset from the training set
Xtr(21:23,:)=[];
Ytr(21:23)=[];

%% Obtain PLS model
[X2,xmean,xstd]=zscore(Xtr);
[Y2,ymean,ystd]=zscore(Ytr); 

Xb=(Xbad-xmean)./xstd;
Xg=(Xgood-xmean)./xstd;

cv=plscv(X2,Y2,5,10);  %10 fold cross validation
lv=cv.optLV;    %Obtain the best LV

[B,W,T,P,Q,Wb,R2X,R2Y]=plsnipals(X2,Y2,lv);

%% Calculate and plot post batch-chart T2
sa=diag(var(T)); %standard deviation
ev=diag(T/sa*T'); %hotelling's values of the training set
CIu = max(bootci(2000,{@median,ev}, 'type', 'per','alpha',.05)); %confidence limits 
%obtain scores
tbad=diag((Xb*W)/sa*(Xb*W)');
tgood=diag((Xg*W)/sa*(Xg*W)');
%plot chart
x = [0 6]; %coordinates for confidence limit
y = [CIu CIu];%coordinates for confidence limit
figure();
subplot(1,2,1)
bar([tbad;tgood])  
title('(a)')
hold on;
line(x,y,'Color','red','LineStyle','--','linewidth',1.5);
ylabel('Hotellings T^2');
xlabel('Batch');
xticks([1:5]);
grid on
box on
xticklabels({'B1','B2','G1','G2','G3'});
legend('T^2','95% Limit');
set(findall(gcf,'-property','FontSize'),'FontSize',24)

%% Calculate and plot post batch-chart SPE
E=X2-T*P';       %Calculate Errors
ev=(diag(E*E'));    %Sum of Errors
CIe = max(bootci(2000,{@median,ev}, 'type', 'per','alpha',.05));%confidence limit
%obtain errors
ebad=(sum((Xb-Xb*W*P').^2,2));
egood=(sum((Xg-Xg*W*P').^2,2));
 %plot chart
x = [0 6]; %coordinates for confidence limit
y = [CIe CIe];%coordinates for confidence limit
subplot(1,2,2)
bar([ebad;egood])
 title('(b)')
hold on;
line(x,y,'Color','red','LineStyle','--','linewidth',1.5);
ylabel('SPE');
xticks([1:5]);
xticklabels({'B1','B2','G1','G2','G3'});
grid on
xlabel('Batch');
legend('SPE','95% Limit');
set(findall(gcf,'-property','FontSize'),'FontSize',24)

%% Calculate and plot normalized cumulative T2
qlim=[];    %Initialize confidence limits
dg1=[];     %Initialize hotelling's value G1
dg2=[];
dg3=[];
db1=[];     %Initialize hotelling's value B1
db2=[];
 for i=1:3:size(X2,2)-1 %Calculate for each hour
    tk =X2(:,i:i+2)*W(i:i+2,:);    %score of training set
    ev=diag(tk/sa*tk');    %Hotelling's of the training set
    CIt = max(bootci(2000,{@median,ev}, 'type', 'per','alpha',.05)) ;%confidence limits
    qlim=[qlim ; CIt];
    %calculate hotelling's values
    tg1 =Xg(1,i:i+2)*W(i:i+2,:);
    dg1=[dg1;  diag(tg1/sa*tg1')];  
    tg2 =Xg(2,i:i+2)*W(i:i+2,:);
    dg2=[dg2;  diag(tg2/sa*tg2')];  
    tg3 =Xg(3,i:i+2)*W(i:i+2,:);
    dg3=[dg3;  diag(tg3/sa*tg3')];   
    tb1 =Xb(1,i:i+2)*W(i:i+2,:);
    db1=[db1; diag(tb1/sa*tb1')];  
    tb2 =Xb(2,i:i+2)*W(i:i+2,:);
    db2=[db2;  diag(tb2/sa*tb2')]; 
 end
%plot chart
figure()
plot(qlim./qlim,'--red','linewidth',1.5);
hold on;
plot ([db1./qlim db2./qlim],'linewidth',1.5);
plot ([dg1./qlim dg2./qlim dg3./qlim],'linewidth',1.5);
legend('95% Limit','B1','B2','G1','G2','G3');
grid on
xlabel('time (h)')
ylabel('Normalized cumulative T^2(k)')
set(findall(gcf,'-property','FontSize'),'FontSize',28)
%% Calculate and plot normalized cumulative SPE
qlim=[];    %Initialize confidence limits
dg1=[];     %Initialize SPE value G1
dg2=[];
dg3=[];
db1=[];     %Initialize SPE value B1
db2=[];
for i=1:3:size(X2,2)
    tk =X2(:,1:i)*W(1:i,:); %calculate scores
    E=X2(:,1:i)-tk*P(1:i,:)'; %Calculate Errors
    ev=diag(E*E'); % SPE
    CIe = max(bootci(2000,{@median,ev}, 'type', 'per','alpha',.05));%confidence limits
    qlim=[qlim ; CIe]; 
    %calculate SPE values
    tg1 =Xg(1,1:i)*W(1:i,:);
    eg1=Xg(1,1:i)-tg1*P(1:i,:)';
    dg1=[dg1;  diag(eg1*eg1')];  
    tg2 =Xg(2,1:i)*W(1:i,:);
    eg2=Xg(2,1:i)-tg2*P(1:i,:)';
    dg2=[dg2;  diag(eg2*eg2')];  
    tg3 =Xg(3,1:i)*W(1:i,:);
    eg3=Xg(3,1:i)-tg3*P(1:i,:)';
    dg3=[dg3;  diag(eg3*eg3')];  
    tb1 =Xb(1,1:i)*W(1:i,:);
    eb1=Xb(1,1:i)-tb1*P(1:i,:)';
    db1=[db1; diag(eb1*eb1')];  
    tb2 =Xb(2,1:i)*W(1:i,:);
    eb2=Xb(2,1:i)-tb2*P(1:i,:)';
    db2=[db2;  diag(eb2*eb2')];     
end
%plot chart
figure ()
plot(qlim./qlim,'--red','linewidth',1.5);
hold on;
plot ([db1./qlim db2./qlim],'linewidth',1.5);
plot ([dg1./qlim dg2./qlim dg3./qlim],'linewidth',1.5);
legend('95% Limit','B1','B2','G1','G2','G3');
grid on
xlabel('time')
ylabel('Normalised cumulative SPE(k)')
set(findall(gcf,'-property','FontSize'),'FontSize',28)

%% Calculate and plot normalized  T2
qlim=[];    %Initialize confidence limits
dg1=[];     %Initialize hotelling's value G1
dg2=[];
dg3=[];
db1=[];     %Initialize hotelling's value B1
db2=[];

 for i=1:1:size(X2,2) %Calculate for each point
    tk =X2(:,i)*W(i,:);    %score of training set
    ev=diag(tk/sa*tk');    %Hotelling's of the training set
    CIt = max(bootci(2000,{@median,ev}, 'type', 'per','alpha',.05)) ;%confidence limits
    qlim=[qlim ; CIt];
    %calculate hotelling's values
    tg1 =Xg(1,i)*W(i,:);
    dg1=[dg1;  diag(tg1/sa*tg1')];  
    tg2 =Xg(2,i)*W(i,:);
    dg2=[dg2;  diag(tg2/sa*tg2')];  
    tg3 =Xg(3,i)*W(i,:);
    dg3=[dg3;  diag(tg3/sa*tg3')];  
    tb1 =Xb(1,i)*W(i,:);
    db1=[db1; diag(tb1/sa*tb1')];  
    tb2 =Xb(2,i)*W(i,:);
    db2=[db2;  diag(tb2/sa*tb2')]; 
 end
%plot chart
figure()
plot(qlim./qlim,'--red','linewidth',1.5);
hold on;
plot ([db1./qlim db2./qlim],'linewidth',1.5);
plot ([dg1./qlim dg2./qlim dg3./qlim],'linewidth',1.5);
legend('95% Limit','B1','B2','G1','G2','G3');
grid on
xlabel('\bf{\it x} (variable x time)')
ylabel('Normalized T^2({\it jk})')
set(findall(gcf,'-property','FontSize'),'FontSize',28)
%% Calculate and plot normalized  SPE
qlim=[];    %Initialize confidence limits
dg1=[];     %Initialize SPE value G1
dg2=[];
dg3=[];
db1=[];     %Initialize SPE value B1
db2=[];
for i=1:1:size(X2,2) %Calculate for each point
    tk =X2(:,i)*W(i,:); %calculate scores
    E=X2(:,i)-tk*P(i,:)'; %Calculate Errors
    ev=diag(E*E'); % SPE
    CIe = max(bootci(2000,{@median,ev}, 'type', 'per','alpha',.05));%confidence limits
    qlim=[qlim ; CIe]; 
    %calculate SPE values
    tg1 =Xg(1,i)*W(i,:);
    eg1=Xg(1,i)-tg1*P(i,:)';
    dg1=[dg1;  diag(eg1*eg1')];  
    tg2 =Xg(2,i)*W(i,:);
    eg2=Xg(2,i)-tg2*P(i,:)';
    dg2=[dg2;  diag(eg2*eg2')];  
    tg3 =Xg(3,i)*W(i,:);
    eg3=Xg(3,i)-tg3*P(i,:)';
    dg3=[dg3;  diag(eg3*eg3')];  
    tb1 =Xb(1,i)*W(i,:);
    eb1=Xb(1,i)-tb1*P(i,:)';
    db1=[db1; diag(eb1*eb1')];  
    tb2 =Xb(2,i)*W(i,:);
    eb2=Xb(2,i)-tb2*P(i,:)';
    db2=[db2;  diag(eb2*eb2')];     
end
%plot chart
figure ()
plot(qlim./qlim,'--red','linewidth',1.5);
hold on;
plot ([db1./qlim db2./qlim],'linewidth',1.5);
plot ([dg1./qlim dg2./qlim dg3./qlim],'linewidth',1.5);
legend('95% Limit','B1','B2','G1','G2','G3');
grid on
xlabel('\bf{\it x} (variable x time)')
ylabel('Normalised SPE({\it jk})')
set(findall(gcf,'-property','FontSize'),'FontSize',28)

