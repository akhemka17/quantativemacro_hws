function vpp1 = evalvp1(cons,x,t,yc)

global betta tetta r g gridx vpfun epsi netw gridy pens pi ny nx 
vpp1 = zeros(ny,1);
for ec=1:ny,
    xp1 = (x-cons)*(1+r)+epsi(t)*netw*gridy(ec)+(1-epsi(t))*pens;
    vpp1(ec) = (func_intp1(gridx,vpfun,xp1));
end;
vpp1 = sum(vpp1'.*pi(yc,:));

end