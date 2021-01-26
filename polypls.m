function [p,q,w,t,u,b,ssqd] = polypls(x,y,lv,n)
% Function to do PLS regression with polynomial inner-relation.
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE
%   p: x-block loading
%   q: y-block loadings
%   w: x-block weights
%   t: x-block scores
%   u: the y-block scores
%   b: matrix of inner-relation coefficients
%   ssqd:variance explained
%   x: The inputs are the matrix of predictor variables
%   y: The vector or matrix of the predicted variable
%   lv: Maximum number of Latent Variables
%   n: Order of the polynomial for the inner-relation
%%
Y=[];
Xs=[];

[mx,nx] = size(x);
[my,ny] = size(y);
p = zeros(nx,lv);
q = zeros(ny,lv);
w = zeros(nx,lv);
t = zeros(mx,lv);
u = zeros(my,lv);
b = zeros(n+1,lv);
ssq = zeros(lv,2);
ssqx = sum(sum(x.^2)');
ssqy = sum(sum(y.^2)');
for i = 1:lv
  [pp,qq,ww,tt,uu] = plsnip(x,y);
  b(:,i) = (polyfit(tt,uu,n))';
  x = x - tt*pp';
  y = y - (polyval(b(:,i),tt))*qq';
  ssq(i,1) = (sum(sum(x.^2)'))*100/ssqx;
  ssq(i,2) = (sum(sum(y.^2)'))*100/ssqy;
  t(:,i) = tt(:,1);
  u(:,i) = uu(:,1);
  p(:,i) = pp(:,1);
  w(:,i) = ww(:,1);
  q(:,i) = qq(:,1);
end
ssqd = zeros(lv,2);
ssqd(1,1) = 100 - ssq(1,1);
ssqd(1,2) = 100 - ssq(1,2);
for i = 2:lv
  for j = 1:2
    ssqd(i,j) = -ssq(i,j) + ssq(i-1,j);
  end
end
