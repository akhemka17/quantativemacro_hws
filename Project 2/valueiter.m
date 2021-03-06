% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% project 2
% solution of household problem for T = \infty


clear all
close all


% -------------------------------------------------------------------------
% SETTINGS
maxit = 50; 
tol=1e-4;
nt=1100;    % periods for simulation
dt=100;     % periods discarded
min_cons=1.0e-08;
crit = 1; % convergence criterion

% parameters
r = 0.02;
rho = 0.03;
g = 0.01;
tetta = 1.2;
betta = 1/(1+rho);
psi=6;
crit = 0.5; 

% grid
nx=20;              % # of grid-points
curv=3.0;           % curvature of grid
xmax = 30;          % scaling factor of saving grid
xmin = -(psi/1+r);
gridx=makegrid(xmin,xmax,nx,curv);
gridx=gridx';

% income shocks
ne = 7;
varepsi = 0.01;
muepsi = -varepsi/2;
% [epsi,probepsi] = qnwnorm(ne,muepsi,varepsi);
% mat=[[1:ne]',epsi,probepsi];
% save rn.txt mat -ascii -double -tabs
mat=load('rn.txt');
epsi=mat(:,2);
probepsi=mat(:,3);
probepsi = repmat(probepsi,1,7);
epsi=exp(epsi);
if (abs(sum(epsi.*probepsi)-1.0)>sqrt(eps)),
    error('random numbers fucked up');
end;
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SOLUTION 
% initial guess


% Discretization of the state space
x = makegrid(xmin,xmax,nx,3)';
xgrid = makegrid(xmin,xmax,nx,3)';
c = zeros(nx,ne);
util1 = zeros(nx,ne);
v = zeros(nx,ne);
dr = zeros(nx,ne);
Tv = zeros(nx,ne);
x = makegrid(xmin,xmax,nx,3)';
gamma=0.1;
iter=0;
while crit>gamma;
    for i=1:nx
    for j=1:ne;
        c = -(((xgrid(i)-epsi(j))/(1+r))*(1+g))+xgrid;
        neg = find(c<=0);
        c(neg) = NaN;
        util1(:,j) = (c.^(1-tetta))/(1-tetta);
        util1(neg,j) = 1.0e-08;
        end
    [Tv(i,:),dr(i,:)] = max(util1+betta*((1+g)^(1-tetta))*(v*probepsi));
    
    end;
    crit = max(max(abs(Tv-v)));
    v = Tv;
    iter = iter+1;
    disp(['no. of iterations', num2str(iter)]);
    end;
    xp = xgrid(dr);
    for j=1:ne;
    c(:,j) = -(((xgrid-epsi(j))/(1+r))*(1+g))+xp(:,j);
    
    end
    
figure;
plot(gridx,v,'LineWidth',2);
ylabel('v');
title('value function');








