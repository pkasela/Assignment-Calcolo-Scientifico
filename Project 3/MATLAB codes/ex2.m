clear all, close all


x0=0;   xL=1;
x=linspace(x0,xL);

mu=1; %coeff di diffusione

%n=[10 20 50 100 200 500 1000];
%bv=[0.1 10 25];

bv=[0.001 0.005 0.01 0.1 1 5 10 50 100 200 300 500];    
n=[200];
u0=1;   uL=0;
Pe=zeros(1,numel(bv));
Pes=zeros(1,numel(bv));
set(0,'DefaultTextInterpreter', 'latex');
for i=1:numel(bv)
    
   for k=1:numel(n)
       clear uh A
       h=abs(xL-x0)/n(k);
       hv(k)=h;
       nnodes=n(k)+1;
       b=bv(i);
       Pe(i)=h*b/(2*mu);
       uex=@(x)(exp(b/mu)-exp(b*x/mu))/(exp(b/mu)-1);

       %muh=mu;               %differenze finite centrate

       % TECNICA UPWIND
       %muh=mu*(1+Pe(i));    %viscosita artificiale
       
       %TECNICA SG
       muh=mu*(Pe(i)+bern(Pe(i)*2)); %viscosita artificiale
       %muh=mu*(Pe(i)+2*Pe(i)/(exp(2*Pe(i))-1));
       
       Pes(i)=h*b/(2*muh);

       %subplot(2,2,i)
       %plot(x,uex(x),'b-'), hold on
       f=zeros(nnodes-2,1);  %nodi interni
       A=zeros(nnodes-2);

       a1=-muh/h^2+b/(2*h);  %colonna j+1
       b1=2*muh/h^2;         %colonna j
       c1=-muh/h^2-b/(2*h);  %colonna j-1

       for j=2:nnodes-3
           A(j,j+1)=a1;
           A(j,j)=b1;
           A(j,j-1)=c1;
       end

       A(1,1)=b1;   A(1,2)=a1;
       A(end,end)=b1;   A(end,end-1)=c1;

       %contributo delle BCs sul termine noto
       f(1)=f(1)-c1*u0;
       f(end)=f(end)-a1*uL;

       uh=A\f;
       uh=[u0;uh;uL];

       xh=[x0:h:xL];
%        plot(xh,uh,'r*-.')
%        str_title=strcat('Pe=',sprintf('%g',Pe(i)),', Pe*=',sprintf('%g',Pes(i)));
%        title(str_title,'fontsize',16)
%        legend('Exact Sol','Numerical Sol','location','SouthWest')
%        xlabel('x','fontsize',14)
%        ylabel('$\varphi$(x)','fontsize',16)
%        hold off
        err(i,k)=max(abs(uh'-uex(xh)));
%        axis ([0 1 0 1.01])
   end
end
%print -dpng Ex2SG_b_25
set(0,'DefaultTextInterpreter', 'tex');
loglog(bv,err(:,1),'*-'),hold on
% loglog(hv,err(2,:),'*-')
% loglog(hv,err(3,:),'*-')
xlabel('\beta')
ylabel('error')
% str_title1=strcat('Error \beta = ',sprintf('%g',bv(1)));
% str_title2=strcat('Error \beta = ',sprintf('%g',bv(2)));
% str_title3=strcat('Error \beta = ',sprintf('%g',bv(3)));
% legend(str_title1,str_title2,str_title3)
str_titlelog=strcat('Error for \Delta x = ',sprintf('%g',hv(1)));
title(str_titlelog)
 %print -dpng Error_SG_n_x_50