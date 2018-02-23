clear 
close

% --- state problem

alpha = 2*10^(-3); %corresponds to diffusion constant

% --- Define constaints and initial condition
L = 1; % length of domain in x direction
tmax = 10; % end time
nx = 100; % number of nodes in x direction
nt = 100; % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);
r = alpha*dt/dx^2; 
r2 = 1 - 2*r;
% --- Loop over time steps
t = 0;
X = 0:dx:L;
%u = normpdf(X,0.3, 10^(-3)); % initial condition
u = sin(pi*X ./ L); % initial condition
u=u';
plot(X,u);
axis([0 1 0 4.5]);
pause(0.001)

% matrix A NOT TIME DEPENDENT
    as =((- alpha)/(dx^2))*ones(1,nx);
    bs = ((1/dt) + (2*alpha/(dx^2))) * ones(1,nx);
    cs = as;

    Diag = [as; bs; cs]';

    d = [-1;0;1];

    A = spdiags(Diag,d,nx,nx);
    
    A(1,1) = 1; A(1,2) =0; % dirichlet boundary condition
    A(end,end) = 1; A(end,end-1) =0; % dirichlet boundary condition

for m=1:nt
uold = u; % prepare for next step
t = t + dt;

    

    % vector b
    b = (1/dt) * uold;
    b(1) = uold(1); b(end) = uold(end); % dirichlet boundary condition

    % calculate system
    
    u=A\b;

plot(X,u);
axis([0 1 0 4.5]);
pause(0.00001)
end