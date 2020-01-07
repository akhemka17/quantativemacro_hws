% COMPLEX VARIANT OF THE KRUSELL SMITH ALGORITHM

function generaleq
opt_det=false;          % 1=deterministic model
opt_nosr=false;         % 1=no survival risk
opt_ny = 2;             % 1=Markov chain with number of states, ny=5,
                        % 2=Markov chain with ny=2 (Krüger-Ludwig calibration)

  
global maxit tol df r nj ny replrate L R delta gridx tau kgrid nk alpha z

tol = 1e-4;     % tolerance level
maxit = 100;    % maximum number of iterations on r
df = 0.1;   

% Calibration:
func_calibr(opt_det,opt_nosr,opt_ny);

% inital guess for psi
psi0 = [0.3;0.3];
psi1 = [1;1];

% approximation for the log of capital
lnk1=zeros(2,nk);
for j=1:nk
    for i=1:2
        lnk1(i,j)=psi0(i)+psi1(i)*log(kgrid(j));
    end
end
k1=exp(lnk1);
tau = func_pens(L,R,replrate);

% solution of household model
[gridx,gridsav,gridasset,cfun,vfun] = func_hh(k1);
end

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function func_calibr(opt_det,opt_nosr,opt_ny)

global betta alpha tetta delta df r nj jr nx nk ny pi gridy pens sr epsi curv pini frac pop totpop grdfac maxit tol L  R replrate z piz kgrid

close all

rho = 0.04;
betta = 1/(1+rho);
tetta = 2;
delta = 0.05;   % depreciation rate of capital
alpha = 0.33;
r = 0.02;

nj=80;
jr=45;
replrate = 0.0;     % pension system replacement rate

nx=30;          % # of grid-points
curv=3.0;       % curvature of grid
grdfac=60;      
   

% deterministic income component:
epsi=ones(nj,1);    
if (jr<nj)
    epsi(jr+1:nj)=0.0;
end

% survival rates
if opt_nosr
    sr = ones(nj,1);
else
    mr = readfile([],'MR.txt',3);
    sr = 1.0-mr(21:21+nj-1,1);
end

% population and fraction living in year...
pop=zeros(nj,1);
pop(1)=100;
for jc=2:nj,
    pop(jc)=pop(jc-1)*sr(jc-1);
end
totpop=sum(pop);

% normalize population to one:
pop=pop/totpop;
totpop=1.0;
frac=pop./totpop;

% working population
L=sum(pop(1:45));
% Aggregate technology shock with two states:
z = [1.03 0.97]';
piz = [0.95 0.05; 0.05 0.95];   % transition matrix

% Grid on aggregate capital: 
Kss = 8;   
nk = 5;
kgrid = linspace(0.5*Kss,1.5*Kss,nk); % end function func_calibr
  
% # of income states
if (opt_det==1)
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
        for tc=1:100
            pini = pi'*pini;
        end
        
        % take exponent and rescale income shocks such that mean is one:
        gridy=exp(gridy)';
        gridy = gridy / sum(gridy.*pini);
        
    else  % markov chain with two states - this is the one we want!
        
        % Alternative -- taken from Krüger and Ludwig (2007) using only two
        % states
        ny = 2;
        
        % transition probability and variance
        rhoeta=0.97;
        vary=0.08;    % taken from Storesletten, Telmer, Yaron
        
        % shock
        epsil=sqrt(vary/(4.0*rhoeta*(1.0-rhoeta)));
        
        % Markov chain
        [pini,pi,gridy]=mchain(rhoeta,epsil);
    end
    
  
end
end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [gridx,gridsav,gridasset,cfun,vfun] = func_hh(k1)

disp('solution of household model');

global alpha betta tetta delta tau nj nx nk ny pi gridy sr epsi curv grdfac replrate kgrid piz z

% grids and decisions rules:
gridx = zeros(nj,ny,nk,2,nx);
gridsav = zeros(nx,1);
gridasset = zeros(nj,ny,nk,2,nx);
cfun = zeros(nj,ny,nk,2,nx);
vfun = zeros(nj,ny,nk,2,nx);
vpfun = zeros(nx,1);
vptrans = zeros(nj,ny,nk,2,nx);

% savings grid: hold it constant:
maxsav=grdfac;
gridsav(2:nx)=makegrid(0.0,grdfac,nx-1,curv);
gridsav(1)=0.0;

% income states
for i=1:ny
    for kgrid=1:nk
        for j=1:2
            
            % current wages and interest rate:
            wage=(1-alpha)*kgrid(kgrid)^(alpha)*z(j);
            ret=alpha*kgrid(kgrid)^(alpha-1)*z(j)-delta;
            
            % cash-on-hand grid at nj:
            inc = epsi(nj)*wage*(1-tau)*gridy(i)+(1-epsi(nj))*replrate*wage*(1-tau);

            % in case of no pension system, assetume some minimum cash on hand:
            minx=max(inc,sqrt(eps));
            maxx=gridsav(nx)*(1.0+ret)+inc;
            gridx(nj,i,kgrid,j,:)=linspace(minx,maxx,nx);

            % Final period Consumption function, assetet holdings, value function, including derivative
            cfun(nj,i,kgrid,j,:)=gridx(nj,i,kgrid,j,:);
            gridasset(nj,i,kgrid,j,:)=(gridx(nj,i,kgrid,j,:)-inc)/(1+ret);
            vfun(nj,i,kgrid,j,:)=U(cfun(nj,i,kgrid,j,:));
            vpfun(:)=MUc(cfun(nj,i,kgrid,j,:));
            vptrans(nj,i,kgrid,j,:)=vpfun.^(-1.0/tetta);
        end
    end
end

% Iterate Backwards
for jc=nj-1:-1:1
    for i=1:ny
        for kgrid=1:nk
            for j=1:2
                kt1=k1(j,kgrid);
                
                %current wage and interest rate
                wage=(1-alpha)*kgrid(kgrid)^(alpha)*z(j);
                ret=alpha*kgrid(kgrid)^(alpha-1)*z(j)-delta;
                
                % interpolate K1 on the grid of aggregate capital
                A = repmat(kt1,[1 length(kgrid)]);
                [xx, stk]=min(abs(A-kgrid));
                K1=kgrid(stk);
                
                for xc=2:nx,
                    vp=zeros(2,2);                
                    for ic=1:ny
                        for jc=1:2
                            % income tomorrow:
                            w1=(1-alpha)*K1^alpha*z(jc);
                            ret1=alpha*K1^(alpha-1)*z(jc)-delta;
                            incp1=epsi(jc+1)*w1*(1-tau)*gridy(ic)+(1-epsi(jc+1))*replrate*w1*(1-tau);

                            % Maximum cash on hand tomorrow:
                            % in case of zero savings and no pension system assetume some
                            % minimum cash on hand
                            cah=max(sqrt(eps),incp1+(1.0+ret1)*gridsav(xc));

                            % Interpolate derivative of value function
                            if ( cah<gridx(jc+1,ic,stk,jc,1))
                                disp('how can this be?')
                            end
                            if ( cah>gridx(jc+1,ic,stk,jc,nx) )
                                % if out of bounds simply set it to decision at nx:
                                vptr = vptrans(jc+1,ic,stk,jc,nx);
                            else
                                vptr = interp1(squeeze(gridx(jc+1,ic,stk,jc,:)),squeeze(vptrans(jc+1,ic,stk,jc,:)),cah);
                            end
                            vp(ic,jc)=vptr.^(-tetta); %value function
                        end %for jc
                    end %for ic
                    
                PI = pi(i,:)'*piz(j,:);
                V = vp.*PI;

                % Euler equation: RHS
                expvp=betta*sr(jc)*(1.0+ret)*sum(sum(V));

                % consumption
                cfun(jc,i,kgrid,j,xc)=invut(expvp);

                % endogenous x-grid:
                gridx(jc,i,kgrid,j,xc)=gridsav(xc)+cfun(jc,i,kgrid,j,xc);
                end %for xc
                    
            % income (wages and pensions) in current period/age:
            inc=epsi(jc)*wage*(1-tau)*gridy(i)+(1-epsi(jc))*replrate*wage*(1-tau);

            % decision at minx
            % notice: correction required for welfare calculation
            % the above is actually slightly inefficient because xmin
            % can be explicitly computed, then gridsav would be age and
            % state dependent.
            minx=max(inc,sqrt(eps));
            if (minx<gridx(jc,i,kgrid,j,2))
                gridx(jc,i,kgrid,j,1)=minx;
            else    % set it to some arbitrary fracion of x(2)
                gridx(jc,i,kgrid,j,1)=0.9*gridx(jc,i,kgrid,j,2);
            end

            % Compute optimal consumption and leisure for minx
            cfun(jc,i,kgrid,j,1)=gridx(jc,i,kgrid,j,1);

            % assetets at all xc:
            gridasset(jc,i,kgrid,j,:)=(gridx(jc,i,kgrid,j,:)-inc)/(1+ret);

            % Update vfun and vpfun
            vpfun(:)=MUc(cfun(jc,i,kgrid,j,:));
            vptrans(jc,i,kgrid,j,:)=vpfun(:).^(-1.0/tetta);

            % Calculate value function
            for xc=1:nx

                v=zeros(2,2);
                for ic=1:ny
                    for jc=1:2
                        % income tomorrow:
                        w1=(1-alpha)*K1^alpha*z(jc);
                        ret1=alpha*K1^(alpha-1)*z(jc)-delta;
                        incp1=epsi(jc+1)*w1*gridy(ic)*(1-tau)+(1-epsi(jc+1))*w1*replrate*(1-tau);

                        % cah tomorrow
                        cah=max(sqrt(eps),incp1+(1.0+ret1)*gridsav(xc));

                        % this should never be the case:
                        if ((cah+0.0001)<gridx(jc+1,ic,stk,jc,1))
                            warning('How can this be ?');
                        end
                        % linear interpolation:
                        v(ic,jc)=func_intp(squeeze(gridx(jc+1,ic,stk,jc,:)),squeeze(vfun(jc+1,ic,stk,jc,:)),cah);
                    end % for jc
                end    % end for ic

                PI=pi(i,:)'*piz(j,:);
                V=v.*PI;
 
                % update value function
                expv=sum(sum(V));
                vfun(jc,i,kgrid,j,xc)=U(cfun(jc,i,kgrid,j,xc))+betta*sr(jc)*expv;
                
                end    % end for xc
            end %for j
        end % for kgrid
    end %for i
end    % end for jc


% ---------------------------------------------------------------------
    function fv = func_intp(x,func,xp)
        
        
        n = length(x);
        if ( xp>x(n) )
            % fv = func(n);
            fv=func_extrapol(x(n-1),x(n),func(n-1),func(n),xp);
        elseif (xp<x(1))
            % fv = func(1);
            fv=func_extrapol(x(1),x(2),func(1),func(2),xp);
        else
            fv = interp1(x,func,xp);
        end
        
    end
% ---------------------------------------------------------------------


% ---------------------------------------------------------------------
    function y=func_extrapol(x1,x2,y1,y2,x)
        
        % simple linear extrapolation
        
        m = (y2-y1)/(x2-x1);
        y = y1 + m*(x-x1);
        
    end
% ---------------------------------------------------------------------

end     % end function func_hh
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [Phi,Phiasset,asset]=func_aggr(gridx,gridsav,cfun,gridasset,wage,ret,tau)
% aggregation and cross-sectional measure
global nj nx ny pi gridy pens sr epsi pini frac totpop replrate

% Compute Cross sectional distributions and aggregate variables
Phi = zeros(nj,ny,nx);          % distribution of assetets conditional by age and shock
Phiasset = zeros(nx,1);             % distribution of assetets

% Distribution of newborns over cash at hand
for i=1:ny
    
    % income (wages and pensions) in current period/age:
    inc=epsi(1)*wage*(1-tau)*gridy(i)+(1-epsi(1))*replrate*wage*(1-tau);
    
    % initial cash-on-hand:
    cahini=inc;
    
    [vals,inds]=basefun(gridx(1,i,:),cahini,nx);
    Phi(1,i,inds(1))=vals(1)*pini(i)*frac(1);
    Phi(1,i,inds(2))=vals(2)*pini(i)*frac(1);
end

for jc=2:nj
    TT = zeros(ny,nx,ny,nx);    % transfer function
    
    for xc=1:nx
        for i=1:ny
            for ic=1:ny
                
                % income (wages and pensions) in current period/age:
                inc=epsi(jc)*wage*(1-tau)*gridy(ic)+(1-epsi(jc))*replrate*wage*(1-tau);
                
                % cash on hand: x=a*(1+r)+y = s(-1)*(1+r)+y;
                cah=inc+(1.0+ret)*gridsav(xc);
                
                [vals,inds]=basefun(gridx(jc,ic,:),cah,nx);
                
                TT(ic,inds(1),i,xc)=vals(1)*pi(i,ic);
                TT(ic,inds(2),i,xc)=vals(2)*pi(i,ic);
            end    
        end    
    end    
   
    for xc=1:nx
        for i=1:ny
            for xcc=1:nx
                for ic=1:ny
                    % transfer distribution:
                    Phi(jc,ic,xcc)=Phi(jc,ic,xcc)+Phi(jc-1,i,xc)*TT(ic,xcc,i,xc)*sr(jc-1);
                end
            end
        end
    end
    
end    % end for jc

% Check that for each country distribution sums to 1
sumprob=sum(sum(sum(Phi(:,:,:))));
if ( ( sumprob < 0.999 ) || ( sumprob > 1.001) )
    beep; beep; beep;
    warning('distribution fucked up');
end

% Check if Grid is Big enough
sumprob=sum(sum(Phi(:,:,nx)));
if (sumprob > 0.001 )
    beep; beep; beep;
    warning('grid too small -- increase your grid');
    pause
end

asset=0.0;
cons=0.0;
lab=0.0;
re=0.0;

% aggregation
for jc=1:nj
    for i=1:ny
        for xc=1:nx
            Phiasset(xc)=Phiasset(xc)+Phi(jc,i,xc);
            
            % assetet holdings = capital stock in general equilibrium
            asset=asset+totpop*Phi(jc,i,xc)*gridsav(xc);
            
            cons=cons+totpop*Phi(jc,i,xc)*cfun(jc,i,xc);
            
            lab=lab+totpop*Phi(jc,i,xc)*gridy(i)*epsi(jc);
            re=re+totpop*Phi(jc,i,xc)*gridy(i)*(1.0-epsi(jc));
        end
    end
end


% ---------------------------------------------------------------------
    function [vals,inds]=basefun(grid_x,x,nx)
        % this subroutine returns the values and the indices of the two basis
        % functions that are positive on a given x in the grid_x
        
        % MF function to lookup the current position
        i=lookup(grid_x,x,0);
        
        if ( (i+1)>nx)
            inds(1)=nx;
            inds(2)=nx;
            vals(2)=0.0;
            vals(1)=1.0;
        elseif (i==0)
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
        end
        
    end 	% end function basefun
% ---------------------------------------------------------------------




end     % end function func_aggr
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [pini,pi,gridy]=mchain(rhoeta,epsil)

% Transition Probabilities
pi=rhoeta*ones(2,2);
pi(1,2)=1.0-rhoeta;
pi(2,1)=1.0-rhoeta;

% Initial Distribution
pini=0.5*ones(2,1);

gridy=zeros(2,1);
gridy(1)=exp(1.0-epsil);
gridy(2)=exp(1.0+epsil);
gridy=2.0*gridy/(sum(gridy));

end  % end function mchain
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++



% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function u = U(c)
global tetta

% utility
if (abs(tetta-1-0)<sqrt(eps))
    u = log(c);
else
    u = c.^(1.0-tetta)/(1.0-tetta);
end
end     % end function U
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function muc=MUc(c)
global tetta

% maringal utility
muc = c.^(-tetta);
end     % end function MUc
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++

% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function invu=invut(marg)
global tetta

% invert utility for c
invu=marg.^(-1.0/tetta);
end     % end function invut
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function grd = makegrid(x1,x2,n,c)
% makes curved grid according to curvature parameter c
scale=x2-x1;
grd(1)=x1;
grd(n)=x2;
for i=2:n-1
    grd(i)=x1+scale*((i-1.0)/(n-1.0))^c;
end
end     % end function makegrid
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function tau = func_pens(L,R,replrate)

tau = replrate*R ./ (L + replrate * R);     % tax on income

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++


% ++++++++++++++++++++++++++++++++++++++++++++++++++++++
function [mpk,Y] = func_firm(asset, L)

global alpha

Y = asset.^alpha * L.^(1-alpha);
ky = asset./Y;
mpk = alpha * ky.^(-1);

end
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++