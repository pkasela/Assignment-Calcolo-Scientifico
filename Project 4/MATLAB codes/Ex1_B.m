clear all, close all

%soluzione con schemi alle DF esplicite del
%pb lineare di trasporto
%du/dt+x*du/dx=0 +u(x,0)=u0(x)=e^(-x^2) + BCs drichlet esatte


T=5;
x0=-5;  xf=5;
dx=1/200;
dt_cr=dx/max(abs(x0),abs(xf)); %valore critico di dt (condiz CFL)

dt=1*dt_cr;

nT=ceil(T/dt);
nX=(xf-x0)/dx;
x=linspace(x0,xf,nX+1);
t=linspace(0,T,nT+1);


%ICs
u0=@(x) 0.*x + exp(-(x).^2).*(and((x>=x0),(x<=xf)));
uex=@(t,x) u0(x.*exp(-t));

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
                u(i,n+1)=u(i,n)-lambda*x(i)/2*(u(i+1,n)-u(i-1,n))...
                    +lambda*abs(x(i))/2*(u(i+1,n)-2*u(i,n)+u(i-1,n));
                %keyboard
            end
            %BCs
         err(n)=(max(abs(uex(n*dt,x)-u(:,n)')));
         end
    case{'LF'}
         disp('metodo--->Lax-Friedrichs');
         for n=1:nT
            u(1,n+1)=uex((n+1)*dt,x0);
            u(nX+1,n+1)= uex((n+1)*dt,xf);
            for i=2:nX
                u(i,n+1)=1/2*(u(i+1,n)+u(i-1,n))...
                    -lambda*x(i)/2*(u(i+1,n)-u(i-1,n));
                %keyboard
            end
            err(n)=(max(abs(uex(n*dt,x)-u(:,n)')));
         end
end %switch
errore=max(err)

time=(0)*dt;
hold on
%traslo avanti di uno perche parto da 1
plot(x,uex(time,x),'g-')
plot(x,u(:,1),'-.','Color',[1 0 0]);
legend('Exact Solution','LF')
xlabel('x')
ylabel('u(x,t)')

time = 1000*dt;
%traslo avanti di uno perche parto da 1
plot(x,uex(time,x),'g-')
plot(x,u(:,1001),'-.','Color',[0.9 0.1 0])
legend('Exact','LF')
% 
time = 2000*dt;
%traslo avanti di uno perche parto da 1
plot(x,uex(time,x),'g-')
plot(x,u(:,2001),'-.','Color',[0.8 0.2 0])
legend('Exact','LF')

time = 3000*dt;
%traslo avanti di uno perche parto da 1
plot(x,uex(time,x),'g-')
plot(x,u(:,3001),'-.','Color',[0.6 0.4 0])
legend('Exact','LF')

% [Y,X]=meshgrid(t,x);
% mesh(Y,X,uex(Y,X)),
% hold on,
% ylabel('X')
% xlabel('T')
% zlabel('U')
