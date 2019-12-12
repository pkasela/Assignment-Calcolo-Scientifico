%% Es 3
% parameters and functions for the HH model of neuron firing  
%uncomment the desired method
clear all, close all
ti=0; tf=100;
ICs=[-65 0 0.3 0.6];
h=0.05;
y(1,:)=ICs;

n=ceil((tf-ti)/h);
t=zeros(1,n);
t(1)=ti;
for i=1:n-1
   t(i+1)=t(i)+h;
   y(i+1,:)=rk4step(@HH,t(i),y(i,:),h);
   %y(i+1,:)=trapstep(@HH,t(i),y(i,:),h);
end

%[t,y]=ode45(@HH,[ti tf],ICs);
%do not ise ode45 while using I_in = 7 for 50<t<=51 since
%it uses bigger steps and does not see the I_in


%to plot all toghether
% figure(1),
% title('Plot using RK4 with all the variables together')
% yyaxis left
% plot(t,y(:,1))
% ylabel('v [mV]')
% yyaxis right
% hold on
% plot(t,y(:,2),'g--')
% plot(t,y(:,3),'r')
% plot(t,y(:,4),'m')
% ylabel('m, n, p')
% xlabel('t [ms]')
% legend('v','m','n','p')

% figure(2),suptitle('Plot using ode45 method and I_{in}(t)=7 for t>50'),
% subplot(2,2,1),plot(t,y(:,1))
% ylabel('v [mV]')
% xlabel('t [ms]')
% subplot(2,2,2),plot(t,y(:,2))
% ylabel('m')
% xlabel('t [ms]')
% subplot(2,2,3),plot(t,y(:,3))
% ylabel('n')
% xlabel('t [ms]')
% subplot(2,2,4),plot(t,y(:,4))
% ylabel('p')
% xlabel('t [ms]')
print -dpng grafico_sovra 
