
function vpp1 = evalvp(cons,x)

global betta tetta r g gridx vpfun epsi probepsi ne nx 

vpp1 = zeros(ne,1);
lambda1 = zeros(ne,1);
for ec=1:ne
    xp1 = (x-cons)*(1+r)/(1+g)+epsi(ec);
    lamda1=vpfun.^(-1/tetta);
    lambda1(ec) =(func_intp(gridx,lamda1,xp1));
    vpp1 = lambda1.^(-tetta);
end
vpp1 = sum(vpp1.*probepsi);
end