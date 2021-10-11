clc
close all
clear

t = -0.5:.001:4.;
dt = t(2) - t(1);
nu = 1e-2;
tic, dom = [-1 1]; x = chebfun('x',dom);

pdefun = @(t,x,u) -u.*diff(u)-diff(u,2)-nu.*diff(u,4);
bc.left = @(u) [u; diff(u)];
bc.right = @(u) [u; diff(u)];

%% Setup preferences for solving the problem.
opts = pdeset('Eps', 1e-4, 'Ylim', [-40,40]);
u0 = exp(-(10*x).^2);
[t, u] = pde15s(pdefun, t, u0, bc, opts);
X = linspace(-1,1,1024);
D = u(X',501:end);

save('Data','D');
