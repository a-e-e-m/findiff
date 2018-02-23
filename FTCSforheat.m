% --- state problem

alpha = 16*10^(-3); %corresponds to diffusion constant

% --- Define constaints and initial condition
L = 1; % length of domain in x direction
tmax = 100; % end time
nx = 100; % number of nodes in x direction
nt = 10000; % number of time steps
dx = L/(nx-1);
dt = tmax/(nt-1);
r = alpha*dt/dx^2; 
r2 = 1 - 2*r;
% --- Loop over time steps
t = 0;
X = 0:dx:L;
u = normpdf(X,0.3, 10^(-3)); % initial condition
% u = sin(pi*X); % initial condition
plot(X,u);
axis([0 1 0 4.5]);
pause(0.001)

for m=1:nt
uold = u; % prepare for next step
t = t + dt;
for i=2:nx-1
u(i) = r*uold(i-1) + r2*uold(i) + r*uold(i+1);
end
plot(X,u);
axis([0 1 0 4.5]);
pause(0.0000001)
end