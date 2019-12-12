function f=HH(t,y)
C=1; 
g1=120;
g2=36;
g3=0.3;

E0=-65;
E1=50;
E2=-77;
E3= -54.4;

alpha_m=@(v) (2.5-0.1*v)./(exp(2.5-0.1*v)-1);
alpha_n=@(v) (0.1-0.01*v)./(exp(1-0.1*v)-1);
alpha_p=@(v) 0.07*exp(-v/20);

beta_m=@(v)4*exp(-v/18);
beta_n=@(v)1/8*exp(-v/80);%changed from 18 to 80 as told in class
beta_p=@(v)1./(exp(3-0.1*v)+1);
f=zeros(1,4); %f=zeros(4,1); for ode45
v=y(1);
m=y(2);
n=y(3);
p=y(4);
%uncomment the Iin desired.
%Iin=0;

% if(t>50) && (t<=51)
%    Iin=7;
% else
%     Iin=0;
% end

% if(t>50)
%    Iin=7;
% else
%     Iin=0;
% end
f(1)=(1/C)*(Iin-g1*m^3*p*(v-E1)-g2*n^4*(v-E2)-g3*(v-E3));
f(2)=(1-m)*alpha_m(v-E0)-m*beta_m(v-E0);
f(3)=(1-n)*alpha_n(v-E0)-n*beta_n(v-E0);
f(4)=(1-p)*alpha_p(v-E0)-p*beta_p(v-E0);
return