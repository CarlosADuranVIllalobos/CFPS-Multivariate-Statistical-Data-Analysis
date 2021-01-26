function [p,q,w,t,u] = plsnip(x,y)
% Function to calculate PLS regressionfor only one LV.
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE
%   p: x-block loadings
%   q: y-block loadings
%   w: x-block weights
%   t: x-block scores
%   u: the y-block scores
%   x: The inputs are the matrix of predictor variables
%   y: The vector or matrix of the predicted variable


[my,ny] = size(y);
if ny > 1
  ssy = sum(y.^2);
  [ymax,yi] = max(ssy);
  u = y(:,yi);
else
  u = y(:,1);
end
conv = 1;
told = x(:,1);
count = 1.0;
%  Specify the conversion tolerance
while conv > 1e-10
  count = count + 1;
  w = (u'*x)';
  w = (w'/norm(w'))';
  t = x*w;
  if ny == 1
    q = 1;
    break
  end
  q = (t'*y)';
  q = (q'/norm(q'))';
  u = y*q;
  conv = norm(told - t);
  told = t;
  
end
p = (t'*x/(t'*t))';
p_norm=norm(p);
t = t*p_norm;
w = w*p_norm;
p = p/p_norm;
