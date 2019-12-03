function generaleq


close all

global nj ny

tic

opt_det=false;          % 1=deterministic model
opt_nosr=false;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,


tol = 1e-4;
maxit = 100;
df = 0.1;

ret = 0.04; 
% SOLUTION

func_calibr(opt_det,opt_nosr,opt_ny);


for it=1:maxit,
    [fval, ass, Y] = func_olg(ret);
    if (abs(fval)<tol),
        disp('convergence');
        break;
    else
        ret = ret - df * fval;
        disp(['iteration #', num2str(it)]);
        disp(['guess of rate of return: ', num2str(ret)]);
    end;
end;
if (it>=maxit),
    warning('no convergence');
end;
disp(['equilibrium interest rate: ', num2str(ret)]);
disp(['equilibrium capital output ratio: ', num2str(ass/Y)]);


end 

function func_calibr(opt_det,opt_nosr,opt_ny)

global betta tetta nj jr nx ny pi gridy netw pens sr epsi curv pini frac pop totpop grdfac delta alpha L

close all

rho = 0.07;
betta = 1/(1+rho);
tetta = 2;
delta = 0.05;
alpha = 0.33;

nj=80;
jr=45;

nx=50;         % # of grid-points
curv=3.0;       % curvature of grid
grdfac=40;      % scaling factor of saving grid

% deterministic income component:
netw=1.0;
pens=0.4;
epsi=ones(nj,1);
if (jr<nj),
    epsi(jr+1:nj)=0.0;
end;

% survival rates
if opt_nosr,
    sr = ones(nj,1);
else
    mr = readfile([],'MR.txt',3);
    sr = 1.0-mr(21:21+nj-1,1);
end;

% population and fraction living in year...
pop=zeros(nj,1);
pop(1)=100;
for jc=2:nj,
    pop(jc)=pop(jc-1)*sr(jc-1);
end;
totpop=sum(pop);

% normalize population to one:
pop=pop/totpop;
totpop=1.0;
frac=pop./totpop;

% working population
L=sum(pop(1:45));


% # of income states
if (opt_det==1),
    ny = 1;
    pini = 1.0;
    gridy = 1.0;
    pi = 1.0;
else
    
    if (opt_ny==1)
        % number of income shocks
        ny = 5;
        % transition probability
        rhoeta=0.98;
        % variance of "permanent" shock
        % taken from Campbell, Viceira, Ch. 7
        vareta=0.01;
        % Markov chain:
        [pi,gridy] = markovappr(rhoeta,sqrt(vareta),2,ny);
        
        % compute invariant distribution
        pini = 1/ny*ones(ny,1);
        for tc=1:100,
            pini = pi'*pini;
        end;
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else
        
        % Alternative -- taken from Kr�ger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        % transition probability and variance
        rhoeta=0.97;
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % shock
        epsil=sqrt(vary/(4.0*rhoeta*(1.0-rhoeta)));
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
    end;
    
end;


end    
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [fval,ass,Y] = func_olg(ret)

global L  delta gridx 


mpk = ret + delta;

wage = func_firm(mpk);

display(wage)

% solution of household model
[gridx,gridsav,gridass,cfun,vfun] = func_hh(wage,ret);

% aggregation
[Phi,PhiAss,ass] = func_aggr(gridx,gridsav,cfun,gridass,wage,ret);


[mpk,Y] = func_mpk(ass, L);

retnew = mpk - delta; 

fval = ret - retnew;




end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [gridx,gridsav,gridass,cfun,vfun] = func_hh(wage,ret)

global betta tetta nj nx ny pi gridy pens sr epsi curv grdfac

disp('solution of household model');

% grids and decisions rules:
gridx = zeros(nj,ny,nx);
gridsav = zeros(nx,1);
gridass = zeros(nj,ny,nx);
cfun = zeros(nj,ny,nx);
vfun = zeros(nj,ny,nx);
vpfun = zeros(nx,1);
vptrans = zeros(nj,ny,nx);

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% income states
for yc=1:ny
    % cash-on-hand grid at nj:
    inc = epsi(nj)*wage*gridy(yc)+(1-epsi(nj))*pens;
    
    % in case of no pension system, assume some minimum cash on hand:
    minx=max(inc,sqrt(eps));
    maxx=gridsav(nx)*(1.0+ret)+inc;
    gridx(nj,yc,:)=linspace(minx,maxx,nx);
    
    % Final period Consumption function, asset holdings, value function, including derivative
    cfun(nj,yc,:)=gridx(nj,yc,:);
    gridass(nj,yc,:)=(gridx(nj,yc,:)-inc)/(1+ret);
    vfun(nj,yc,:)=U(cfun(nj,yc,:));
    vpfun(:)=MUc(cfun(nj,yc,:));
    vptrans(nj,yc,:)=vpfun.^(-1.0/tetta);
end;

% Iterate Backwards
for jc=nj-1:-1:1,
    
    for yc=1:ny
        
        for xc=2:nx,
            vp=zeros(2,1);
            
            for ycc=1:ny,
                % income tomorrow:
                incp1=epsi(jc+1)*wage*gridy(ycc)+(1-epsi(jc+1))*pens;
                
                % Maximum cash on hand tomorrow:
                % in case of zero savings and no pension system assume some
                % minimum cash on hand
                cah=max(sqrt(eps),incp1+(1.0+ret)*gridsav(xc));
                
                % Interpolate derivative of value function
                if ( cah<gridx(jc+1,ycc,1)),
                    disp('how can this be?')
                end;
                if ( cah>gridx(jc+1,ycc,nx) ),
                    % if out of bounds simply set it to decision at nx:
                    vptr = vptrans(jc+1,ycc,nx);
                else
                    vptr = interp1(squeeze(gridx(jc+1,ycc,:)),squeeze(vptrans(jc+1,ycc,:)),cah);
                end;
                vp(ycc)=vptr.^(-tetta);
            end;
            
            % Euler equation: RHS
            expvp=betta*sr(jc)*(1.0+ret)*sum(pi(yc,:)*vp(:));
            
            % consumption
            cfun(jc,yc,xc)=invut(expvp);
            
            % endogenous x-grid:
            gridx(jc,yc,xc)=gridsav(xc)+cfun(jc,yc,xc);
        end;
        
        % income (wages and pensions) in current period/age:
        inc=epsi(jc)*wage*gridy(yc)+(1-epsi(jc))*pens;
        
        % decision at minx
        % notice: correction required for welfare calculation
        % the above is actually slightly inefficient because xmin
        % can be explicitly computed, then gridsav would be age and
        % state dependent.
        minx=max(inc,sqrt(eps));
        if (minx<gridx(jc,yc,2)),
            gridx(jc,yc,1)=minx;
        else    % set it to some arbitrary fracion of x(2)
            gridx(jc,yc,1)=0.9*gridx(jc,yc,2);
        end;
        
        % Compute optimal consumption and leisure for minx
        cfun(jc,yc,1)=gridx(jc,yc,1);
        
        % assets at all xc:
        gridass(jc,yc,:)=(gridx(jc,yc,:)-inc)/(1+ret);
        
        % Update vfun and vpfun
        vpfun(:)=MUc(cfun(jc,yc,:));
        vptrans(jc,yc,:)=vpfun(:).^(-1.0/tetta);
        
        % Calculate value function
        for xc=1:nx,
            
            v=zeros(2,1);
            for ycc=1:ny,
                % income tomorrow:
                incp1=epsi(jc+1)*wage*gridy(ycc)+(1-epsi(jc+1))*pens;
                
                % cah tomorrow
                cah=max(sqrt(eps),incp1+(1.0+ret)*gridsav(xc));
                
                % this should never be the case:
                if ((cah+0.0001)<gridx(jc+1,ycc,1)),
                    warning('How can this be ?');
                end;
                % linear interpolation:
                v(ycc)=func_intp(squeeze(gridx(jc+1,ycc,:)),squeeze(vfun(jc+1,ycc,:)),cah);
            end;    % end for ycc
            
            % update value function
            expv=sum(pi(yc,:)*v(:));
            vfun(jc,yc,xc)=U(cfun(jc,yc,xc))+betta*sr(jc)*expv;
        end;    % end for xc
        
    end;    % end for yc
    
end;    % end for jc


end     % end function func_hh
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

function [Phi,PhiAss,ass]=func_aggr(gridx,gridsav,cfun,gridass,wage,ret)

global  nj nx ny pi gridy pens sr epsi pini frac totpop 

disp('aggregation and cross-sectional measure');

% Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx);          % distribution of assets conditional by age and shock
PhiAss = zeros(nx,1);             % distribution of assets

% Distribution of newborns over cash at hand
for yc=1:ny
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*wage*gridy(yc)+(1-epsi(1))*pens;
    
    % initial cash-on-hand:
    cahini=inc;
    
    [vals,inds]=basefun(gridx(1,yc,:),cahini,nx);
    Phi(1,yc,inds(1))=vals(1)*pini(yc)*frac(1);
    Phi(1,yc,inds(2))=vals(2)*pini(yc)*frac(1);
end;

for jc=2:nj
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
        for yc=1:ny
            for ycc=1:ny
                
                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*wage*gridy(ycc)+(1-epsi(jc))*pens;
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+ret)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ycc,:),cah,nx);
                
                TT(ycc,inds(1),yc,xc)=vals(1)*pi(yc,ycc);
                TT(ycc,inds(2),yc,xc)=vals(2)*pi(yc,ycc);
            end;    
        end;    
    end;    
    
    for xc=1:nx
        for yc=1:ny
            for xcc=1:nx
                for ycc=1:ny
                    % transfer distribution:
                    Phi(jc,ycc,xcc)=Phi(jc,ycc,xcc)+Phi(jc-1,yc,xc)*TT(ycc,xcc,yc,xc)*sr(jc-1);
                end;
            end;
        end;
    end;
    
end;    % end for jc


% Check that for each country distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) ),
    beep; beep; beep;
    warning('distribution fucked up');
end;

% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 ),
    beep; beep; beep;
    warning('grid too small -- increase your grid');
    pause
end;

ass=0.0;
cons=0.0;


% aggregation
for jc=1:nj
    for yc=1:ny
        for xc=1:nx,
            PhiAss(xc)=PhiAss(xc)+Phi(jc,yc,xc);
            
            % asset holdings = capital stock in general equilibrium
            ass=ass+totpop*Phi(jc,yc,xc)*gridsav(xc);
            
            cons=cons+totpop*Phi(jc,yc,xc)*cfun(jc,yc,xc);
            
       
%             
        end;
    end;
end;


% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookup1(grid_x,x,0);
        
        if ( (i+1)>nx),
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0),
            inds(1)=1;
            inds(2)=1;
            vals(1)=1.0;
            vals(2)=0.0;
        else
            inds(1)=i;
            inds(2)=i+1;
            dist = grid_x(i+1)-grid_x(i);
            vals(2)=( x-grid_x(i) )/dist;
            vals(1)=( grid_x(i+1)-x )/dist;
        end;
        
    end 	% end function basefun
% ---------------------------------------------------------------------




end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
