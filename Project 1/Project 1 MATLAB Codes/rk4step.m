function y=rk4step(odefun,t,un,h)

K1=odefun(t,un);%passo di EE
K2=odefun(t+h/2,un+h*K1/2);
K3=odefun(t+h/2,un+h*K2/2);
K4=odefun(t+h,un+h*K3);

y=un+h*(K1+2*K2+2*K3+K4)/6;
return