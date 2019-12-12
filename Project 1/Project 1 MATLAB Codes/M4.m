%% metodo M4 restituisce u_{n+1}
function yn1=M4(odefun,tn,yn,h)
% e' un metodo RK a 3 stadi
K1=odefun(tn,yn);
tK2=tn+(1/3)*h;
yK2=yn+(1/3)*h*K1;
K2=odefun(tK2,yK2);
tK3=tn+(2/3)*h;
yK3=yn+(2/3)*h*K2;
K3=odefun(tK3,yK3);
yn1=yn+(h/4)*(K1+K3);

return