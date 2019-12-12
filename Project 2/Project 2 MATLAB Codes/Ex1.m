%% Exercise 1
clear all, close all

xi=0; xf=1;
nn=[15];%thus total nodes are 16
%% Problem f(x)=1

alpha_ex=@(x) x;
beta_ex=@(x) x.^2/2;
f=@(x) 1;
u_ex=@(x) (x-x.^2)/2;

%% Problem f(x)=x

% alpha_ex=@(x) x.^2/2;
% beta_ex=@(x) x.^3/3;
% f=@(x) x;
% u_ex=@(x) (x-x.^3)/6;

%% Problem f(x)=exp(x)

% alpha_ex=@(x) exp(x)-1;
% beta_ex=@(x) exp(x).*(x-1)+1;
% f=@(x) exp(x);
% u_ex=@(x) 1+(exp(1)-1)*x-exp(x);

%% main
alpha(1)=0;
beta(1)=0;
x(1)=0;
u(1)=0;

z=linspace(0,1);
for k=1:numel(nn)
    n=nn(k);
    h=(xf-xi)/n;
    hh(k)=h;
    for i=1:n

        x(i+1)=x(i)+h;
        x12=x(i)+h/2;
        alpha(i+1)=alpha(i) + h/4*(f(x(i))+2*f(x12)+f(x(i+1)));
        beta(i+1)=beta(i)+h/4*(x(i)*f(x(i))+2*x12*f(x12)+x(i+1)*f(x(i+1)));
    end

    for i=1:n
        u(i+1)=x(i+1)*(alpha(end)-beta(end))+beta(i+1)-x(i+1)*alpha(i+1);
    end

    % plot(z,alpha_ex(z),'r'), hold on
    % plot(x,alpha,'g')
    % 
    % figure(2), plot(z,beta_ex(z),'r'), hold on
    % plot(x,beta,'g')

    plot(z,u_ex(z),'r'), hold on
    plot(x,u,'g-.')
    legend('Exact Solution','Numerical Solution')
    %print -dpng Solf_exp(x) %to save the image
    err(k)=max(abs(u_ex(x)-u));
end
close all
loglog(hh,err,'-*')
legend('error for f=e^x')
xlabel('\Delta x')
ylabel('Error')
