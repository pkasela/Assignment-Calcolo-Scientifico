%% Exercise 2
clear all, close all
%% Problem details
xi=0; xf=1;
theta=1/2;
ti=0; tf=1;

epsilon=@(x) 2+cos(pi*x);
uex=@(x,t) x+sin(pi*x).*exp(-t);
uIC=@(x) x+sin(pi*x);
u0=0; uL=1; %\forall t
%NNx=[50 100 200 400 800];
NNy=[10 20 40 80 160 320];
for l=1:numel(NNy)
%Nx=NNx(l); %node on x
Nx=2000;
Ny=NNy(l); %nodes on t
nx=Nx-1; %internal nodes
h=(xf-xi)/Nx; %delta x
%hh(l)=h;
ht=(tf-ti)/Ny; %delta t
hh(l)=ht;
clear u
x=linspace(xi,xf,Nx+1);
t=linspace(ti,tf,Ny+1);
u(:,1)=uIC(x(2:end-1));
%% calculation of f(x,t)
% syms u(x,t) eps(x)
% 
% u=x+sin(pi*x).*exp(-t);
% eps=2+cos(pi*x);
% 
% ut=diff(u,'t');
% ux=diff(u,'x');
% f=ut-diff(eps.*ux,'x');
% f=simplify(f);
f=@(x,t) exp(-t).*sin(pi*x).*(pi*exp(t) + 2*pi^2*cos(pi*x) + 2*pi^2 - 1);
    %% create matrix A (indipendent from time)
    A=zeros(nx);

    for j=2:nx-1
        xm12=(j*h-h/2); %node x_{j-1/2}
        xp12=(j*h+h/2); %node x_{j+1/2}
        A(j,j-1)=-epsilon(xm12);
        A(j,j)=(epsilon(xm12)+epsilon(xp12));
        A(j,j+1)=-epsilon(xp12);
    end

    %first row of matrix A
    xm12=h/2; xp12=3/2*h; 
    A(1,1)=epsilon(xm12)+epsilon(xp12);
    A(1,2)=-epsilon(xp12);

    %last row of matrix A
    xm12=xf-3/2*h; xp12=xf-h/2; 
    A(end,end)=epsilon(xm12)+epsilon(xp12);
    A(end,end-1)=-epsilon(xm12);

    A=A/h^2;

    Atilde=eye(nx)/ht+theta*A;

%% calculate u
for k=2:numel(t)
    %% known term b

    fk=zeros(nx,1);
    for i=1:nx
        fk(i)=f(i*h,t(k-1));
    end

    fk1=zeros(nx,1);
    for i=1:nx
        fk1(i)=f(i*h,t(k));
    end

    u_remaining=zeros(nx,1);
    u_remaining(1)=(1-theta)*u0/h^2*epsilon(xm12)+theta*u0/h^2*epsilon(xm12);
    u_remaining(end)=(1-theta)*uL/h^2*epsilon(xp12)+theta*uL/h^2*epsilon(xp12);

    btilde = (eye(nx)/ht-(1-theta)*A)*u(:,k-1)+theta*fk1+(1-theta)*fk+u_remaining;

    %% calculate u_{k+1}
    u(:,k)=Atilde\btilde;
end
%% plot of u

 u=[u0*ones(1,Ny+1);u;uL*ones(1,Ny+1)];
 x=xi:h:xf; t=ti:ht:tf;
 [T,X]=meshgrid(t,x);
 surf(T,X,uex(X,T))
 title('Exact solution')
 xlabel('t')
 ylabel('x')
 zlabel('u(x,t)')
 %print -dpng Ex2ExSol
 figure, surf(T,X,u)
  title('Numerical solution')
 xlabel('t')
 ylabel('x')
 zlabel('u(x,t)')
 %print -dpng Ex2NumSol


 error(l)=max(max(abs(u-uex(X,T))));
end

figure, loglog(hh,error,'-*')
xlabel('\Deltat')
ylabel('Error')
title('Error for \Deltax fixed')
legend('Error line')