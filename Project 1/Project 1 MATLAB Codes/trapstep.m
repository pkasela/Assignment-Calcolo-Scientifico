function y=trapstep(odefun,t,z,h)
    
%passo predictor(passo EE)
z1=z+h*odefun(t,z);

%passo corrector
z2=odefun(t+h,z1);
%media
y=z+(odefun(t,z)+z2)*h/2;
return