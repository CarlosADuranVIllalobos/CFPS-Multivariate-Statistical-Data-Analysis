
function [B,limit] = clim(X,Y,lv)
% Function to do PLS regression and confidence limits using bootstrap calculation 
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE

%   B: Regression vector
%   limit: Upper confidence limit
%   X: Predictor variables
%   Y: Response variables
%   lv=Number of Latent Variables
xmat=size(X,2);
ymat=size(Y,2);

B = plsnipals(X,Y,lv);
n=50000;

[x,bootsamp] = bootstrp(n,@(bootr)plsnipals(bootr(:,1:xmat),bootr(:,xmat+ymat),lv),[X,Y]);%bootstrap calculation

se=0;    
for i=1:1:n
se=se+(B'-x(i,:)).^2;
end
see1=sqrt(se/(n-1));
see=see1.*2;

limit = see;
end

