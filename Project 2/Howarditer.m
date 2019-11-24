% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% project 2
% solution of household problem for T = \infty

clear all
close all

% -------------------------------------------------------------------------
% SETTINGS
maxit = 100; 
tol=1e-4;
nt=1100;    % periods for simulation
dt=100;     % periods discarded
min_cons=1.0e-08;

% parameters
r = 0.02;
rho = 0.03;
g = 0.01;
tetta = 2;
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
if (abs(sum(epsi.*probepsi)-1.0)>sqrt(epsi))
    error('random numbers fucked up');
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% SOLUTION 
% initial guess

% Discretization of the state space
xgrid = makegrid(xmin,xmax,nx,3)';
%I change also here xpo...since I follow the document
xp0 = makegrid(xmin,xmax,nx,3)';
c = zeros(nx,ne);
util1 = zeros(nx,ne);
v = zeros(nx,ne);
Tv = zeros(nx,ne);
EV=zeros(nx,ne);
dr0 = repmat([1:nx]',1,ne); % initial guess
dr=zeros(nx,ne);
iter=0;
gamma=0.01;

while crit>gamma
    for i=1:nx
    for j=1:ne
        c = -(((xgrid(i)-epsi(j))/(1+r))*(1+g))+xgrid;
        neg = find(c<=0);
        c(neg) = NaN;
        util1(:,j) = (c.^(1-tetta))/(1-tetta);
        util1(neg,j) = 1.0e-08;
        EV(:,j) = v*probepsi(j,:)';
        end
     [v1,dr(i,:)] = max(util1+betta*((1+g)^(1-tetta))*EV);
    end
    % decision rules
    xp=xgrid(dr);
    Q=sparse(nx*ne,nx*ne);
    for j=1:ne
        c= -(((xgrid(i)-epsi(j))/(1+r))*(1+g))+xp(:,j);
      % update the value
        util1(:,j) = (c.^(1-tetta))/(1-tetta);
        Q0=sparse(nx,nx);
        for i=1:nx
            Q0(i,dr(i,j)) = 1;
        end
        Q((j-1)*nx+1:j*nx,:) = kron(probepsi(j,:),Q0);
    end
    Identity=speye(nx*ne);
    TV=(Identity-betta*Q)\util1(:);
    crit = max(max(abs(xp-xp0)));
    v = reshape(TV,nx,ne);
    xp0=xp;
    iter = iter+1;
    disp(['no. of iterations', num2str(iter)]);
    end;
    c=zeros(nx,ne);
for j=1:ne
    c= -(((xgrid(i)-epsi(j))/(1+r))*(1+g))+xp(:,j);
end
figure;
plot(gridx,v,'LineWidth',2);
ylabel('v');
title('value function');


function grd = makegrid(x1,x2,nc,curv)

% makes curved grid according to curvature parameter c
scale=x2-x1;
grd(1)=x1;
grd(nc)=x2;
for ic=2:nc-1
    grd(ic)=x1+scale*((ic-1.0)/(nc-1.0))^curv;
end

end     






