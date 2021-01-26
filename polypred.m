function y = polypred(x,b,p,q,w,lv)
% Function to predict using polynomial PLS.
%% Copyright
% Carlos Alberto Duran-Villalobos June 2020 University of Manchester.
% Data provided by UCL and Sutro
% Copyright (c) Future Targeted Healthcare Manufacturing Hub
% Reference: "Multivariate statistical data analysis of cell-free protein synthesis towards monitoring and control", AIChE
%   y: the predicted variable
%   p: Outputs are the x-block loading
%   q: y-block loadings
%   w: x-block weights
%   t: x-block scores
%   b: matrix of inner-relation coefficients
%   x: The inputs are the matrix of predictor variables
%   lv: Number of Latent Variables

%%

[mx,nx] = size(x);
[mq,nq] = size(q);
[mw,nw] = size(w);
that = zeros(mx,lv);
y = zeros(mx,mq);

%  Start by calculating all the xblock scores
for i = 1:lv
  that(:,i) = x*w(:,i);
  x = x - that(:,i)*p(:,i)';
end
%  Use the xblock scores and the b to build up the prediction
for i = 1:lv
y = y + (polyval(b(:,i),that(:,i)))*q(:,i)';
end
