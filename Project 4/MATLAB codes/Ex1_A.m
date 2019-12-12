clear all, close all

%soluzione con schemi alle DF esplicite del
%pb lineare di trasporto
%du/dt+a*du/dx=0 +u(x,0)=u0(x) + BCs , a=2

a=2; %velocita di propagazione

T=5;
x0=-5;  xf=5;
dx=1/25;
dt_cr=dx/abs(a); %valore critico di dt (condiz CFL)

dt=1*dt_cr;

nT=ceil(T/dt);
nX=(xf-x0)/dx;
x=linspace(x0,xf,nX+1);
t=linspace(0,T,nT+1);


%ICs
u0=@(x) 0.*x + exp(-(x).^2).*(and((x>=x0),(x<=xf)));
uex=@(t,x) u0(x-a*t);

lambda=dt/dx;
u=zeros(nX+1,nT+1);
u(:,1)=u0(x); %ICs
metodo='LF';
switch (metodo)
    case{'upwind'}
        disp('metodo--->upwind');
        for n=1:nT
             u(1,n+1)=uex((n+1)*dt,x0);
             u(nX+1,n+1)=uex((n+1)*dt,xf);
            for i=2:nX
                u(i,n+1)=u(i,n)-lambda*a/2*(u(i+1,n)-u(i-1,n))...
                    +lambda*abs(a)/2*(u(i+1,n)-2*u(i,n)+u(i-1,n));
            end 
         err(n)=(max(abs(uex(n*dt,x)-u(:,n)')));
         end
    case{'LF'}
        disp('metodo--->Lax-Friedrichs');
        for n=1:nT
            u(1,n+1)=uex((n+1)*dt,x0);
            u(nX+1,n+1)= uex((n+1)*dt,xf);
            for i=2:nX
                u(i,n+1)=1/2*(u(i+1,n)+u(i-1,n))...
                    -lambda*a/2*(u(i+1,n)-u(i-1,n));
            end
            err(n)=(max(abs(uex(n*dt,x)-u(:,n)')));
         end
end %switch
errore=max(err)
% plot 
figure,
subplot(2,2,1)
time=(0)*dt;
 plot(x,uex(time,x),'m-'), hold on
 %traslo avanti di uno perche parto da 1
 plot(x,u(:,1),'b-.');
 legend('Exact for t=0','LF')
 xlabel('x')
 ylabel('u(x,t)')
 title('Plot for t=0')
 
 subplot(2,2,2)
time=(50)*dt;
 plot(x,uex(time,x),'m-'), hold on
 %traslo avanti di uno perche parto da 1
 plot(x,u(:,51),'b-.');
 legend('Exact for t=0.5','LF')
 xlabel('x')
 ylabel('u(x,t)')
 title('Plot for t=0.5')
 
subplot(2,2,3)
title('Plot for t=1')
 time=(100)*dt;
 plot(x,uex(time,x),'m-'), hold on
 %traslo avanti di uno perche parto da 1
 plot(x,u(:,101),'b-.');
 legend('Exact for t=1','LF')
 xlabel('x')
 ylabel('u(x,t)')
 title('Plot for t=1')
 
 time=(200)*dt;
 subplot(2,2,4)
 
 plot(x,uex(time,x),'m-'), hold on
 %traslo avanti di uno perche parto da 1
 plot(x,u(:,201),'b-.');
 legend('Exact for t=2','LF')
 xlabel('x')
 ylabel('u(x,t)')
 title('Plot for t=2')

%  [Y,X]=meshgrid(t,x);
%  mesh(Y,X,uex(Y,X)),
%  title('Exact Solution')
%  ylabel('X')
%  xlabel('T')
%  zlabel('U(t,x)')
%  figure,
%  mesh(Y,X,u)
%  title('Numerical Solution')
%  ylabel('X')
%  xlabel('T')
%  zlabel('U(t,x)')
 
% figure, hold on
% plot(x,x/a,'b')
% for i=1:10
%     plot(x,x/a+3*i/2,'r')
%     plot(x,x/a-i*3/2,'g')
% end
% axis([-10 10 0 5])
% xlabel('x')
% ylabel('t')
% title('characterisitic lines for the method (A)')
