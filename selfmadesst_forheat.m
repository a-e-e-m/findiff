clear
close

% --- steady state heat equation in R^2

alpha = 16*10^(-2); %diffusion constant


% --- Define constaints and initial condition
Lx = 1; % length of domain in x direction
Ly = 1; % length of domain in y direction

dxy = 0.001; % space step

Nx = Lx/dxy + 1; % number of nodes in x direction
My = Ly/dxy + 1; % number of nodes in y direction

u = zeros(My,Nx); % initialise solution u
u(ceil(My/2) -10 : ceil(My/2) +10, 1) = 1; % initial condition


% [X1,X2] = meshgrid(linspace(0,Lx,Nx)',linspace(0,Ly,My)');
% u = Conc_Gauss( X1, X2, 0, 0.5 ,0.01, 0.01 );
% % surf(X1, X2, Y)

% matrix acting on X
    dim = Nx*My;

    a = ones(dim,1);
    b = ones(dim,1);
    c = 4 * ones(dim,1);
    d = ones(dim,1);
    e = ones(dim,1);

    Diag = [a,b,c,d,e];

    d = [-My;-1;0;1;My];

    A = spdiags(Diag,d,dim,dim);
    
    for k=-10:1:10;
        s=ceil(My^2 /2 + k*My +1);
        A(s, :) = sparse(1,[s],1,1,dim);
    end
    
% b in AX=b
    b = u(:); % converting matrix to column vector
    b = ((-1)/alpha) * b;
    
% solve system

X = A\b;

% plot
U = reshape(X, [My Nx]);

x=0:dxy:Lx;
y=0:dxy:Ly;
[xx,yy] = meshgrid(x,y);

surf(xx,yy,U,'edgecolor','none')
