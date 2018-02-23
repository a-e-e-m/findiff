% the equation to solve is - (d^2 c(x))/(dx)^2 = R(x)
% where c is a function from 0:1 to R
%       R(x) is 1 for x=0.3 and 0 else
%       c(0)=c(1)=0
% method is forward euler

clear

delta = 0.01; % space step
N = 1/delta; % number of steps
samplepoints = 0:delta:1;

c = zeros(N+1,1);
c(1,1) = 0;
c(end,1) = 0;

% R = zeros(N+1,1);
% peak = 0.3 / delta;
% R(peak) = 1;


R = normpdf(samplepoints,0.3, 10^(-3));
R = mat2gray(R);
plot(samplepoints, R);


% define matrix equation coming from forward euler method

    % matrix A
    Diagminus1 =[ (-0.5)*ones(1,N-2), 0];
    Diag0 = ones(1,N-1);
    Diagplus1 = [0, (-0.5)*ones(1,N-2)];

    Diag = [Diagminus1; Diag0; Diagplus1]';

    d = [-1;0;1];

    A = spdiags(Diag,d,N-1,N-1);

    % vector b
    b = 0.5 * delta^2 * R(1,2:N)';
    
    boundarycond = zeros(N-1,1);
    boundarycond(1,1) = c(1,1);
    boundarycond(end,1) = c(end,1);
    
    b = b + boundarycond;
    
    
% compute
X = A\b;

% solution
c(2:end-1) = X;



plot(samplepoints, c)

