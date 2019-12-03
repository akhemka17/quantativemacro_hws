% Consumption equivalent 

clear, clc

v0 = importdata('vT_dec.txt');  
v1 = importdata('vT1_dec.txt'); 
xT = importdata('xT_dec.txt'); 

theta = 1.1;

% CEV:
g = ((v1./v0).^(1/(1-theta)))-1

plot(xT(2,:),g(2,:))