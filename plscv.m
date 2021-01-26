function CV=plscv(X,Y,lv,K)
% Function to do K-fold Cross-validation for PLS 
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE
%   CV: Output structure
%   Y: Matrix of the predicted variables
%   X: matrix of predictor variables
%   lv: Number of Latent Variables
%   K: fold. 



if nargin<4;K=10;end
if nargin<3;lv=3;end

[Y,indexyy]=sort(Y);
X=X(indexyy,:);

[Mx,Nx]=size(X);
lv=min([size(X) lv]);
ystest=nan(Mx,1);
Yr=nan(Mx,lv);

groups = 1+rem(0:Mx-1,K);
for group=1:K
    calk = find(groups~=group);
    testk = find(groups==group);  
    Xcal=X(calk,:);
    ycal=Y(calk);
    Xtest=X(testk,:);
    ytest=Y(testk);
    
    %  data pretreatment
    [Xs,xmean,xstd]=zscore(Xcal);
    [ys,ymean,ystd]=zscore(ycal);   

    [B,Wstar,T,P,Q]=plsnipals(Xs,ys,lv);   
 
    yp=[];
    for j=1:lv
        B=Wstar(:,1:j)*Q(1:j);
        C=ystd*B./xstd';
        coef=[C;ymean-xmean*C;];
        Xteste=[Xtest ones(size(Xtest,1),1)];
        ypred=Xteste*coef;
        yp=[yp ypred];
    end
    
    Yr(testk,:)=[yp];
    ystest(testk,:)=ytest;
end

%	 return the original order
Yr(indexyy,:)=Yr;
Y(indexyy)=Y;

%	 mean and sd of squared error 
error=Yr-repmat(Y,1,lv);
error2=error.^2;
error2_m=sum(error2)/Mx;
error2_std= sqrt(sum((error2-repmat(mean(error2),Mx,1)).^2)/(Mx-1)); % unbiased estimator

  %	  Q2
cv=sqrt(error2_m);
[RMSEP,index]=min(cv);index=min(index);
SST=sumsqr(ystest-mean(Y));
for i=1:lv
    SSE=sumsqr(Yr(:,i)-Y);
    Q2(i)=1-SSE/SST;
end
  
  indexSD=find(error2_m <= min(error2_m)+error2_std(index));
  indexSD=min(indexSD);

%   output 
  CV.Ypred=Yr;
  CV.predError=error;
  CV.RMSECV=cv;
  CV.Q2=Q2;
  CV.RMSECV_min=RMSEP;
  CV.Q2_max=Q2(index);
  CV.optLV=index;
  CV.RMSECV_min_1SD=cv(indexSD);
  CV.Q2_max_1SD=Q2(indexSD);
  CV.optLV_1SD=indexSD;






