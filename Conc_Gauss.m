function [ y ] = Conc_Gauss( x1, x2, h1, h2,sigma1, sigma2 )
%input is a point (x1,x2) and the souce location (h1,h2) in R^2
%output is the concentration at point (x1,x2) in (0,1) 
%assuming a two-dimensional Gaussian function with center in (h1,h2) and 
%standard deviation 1 as well as factor 1
% see https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function

A=1; %factor
y = A* exp( -( (x1-h1).^2 / (2*sigma1.^2)  + (x2-h2).^2 / (2*sigma2.^2) ) );
end

