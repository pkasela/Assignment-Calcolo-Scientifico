%% Exercise 3
clear all, close all

xi=0; xf=1;
dx=0.05;
x=xi:dx:xf;

ti=0; tf=15;
dt=0.8*0.00125;
t=ti:dt:tf;


%ICs
%u_in= 1*(ones(numel(x),1));
u_in=-x.^3'+(3/2)*x.^2';
%v_in= 0*(ones(numel(x),1));
v_in=cos(pi*x).^2';

%  v_in(1)=2;
%  v_in(end)=2;

u=zeros(numel(t),numel(x));
v=zeros(numel(t),numel(x));

u(1,:)=u_in;
v(1,:)=v_in;
for k=1:numel(t)-1
    for j=2:numel(x)-1
        u(k+1,j)=u(k,j)+dt/dx^2*(u(k,j-1)-2*u(k,j)+u(k,j+1))+dt*(u(k,j)-u(k,j)*v(k,j));
        v(k+1,j)=v(k,j)+dt/dx^2*(v(k,j-1)-2*v(k,j)+v(k,j+1))+dt*(-v(k,j)+u(k,j)*v(k,j));
    end
    u(k+1,1)=u(k,1)+dt/dx^2*(2*u(k,2)-2*u(k,1))+dt*u(k,1)*(1-v(k,1)); 
    u(k+1,end)=u(k,end)+dt/dx^2*(2*u(k,end-1)-2*u(k,end))+dt*u(k,end)*(1-v(k,end));
    v(k+1,1)=v(k,1)+dt/dx^2*(2*v(k,2)-2*v(k,1))+dt*v(k,1)*(u(k,1)-1); 
    v(k+1,end)=v(k,end)+dt/dx^2*(2*u(k,end-1)-2*u(k,end))+dt*v(k,end)*(u(k,end)-1);
end

%% mesh plot
[T,X]=meshgrid(t,x);
figure, subplot(2,1,1)
mesh(T,X,u')
xlabel('time')
ylabel('space')
zlabel('u population')
subplot(2,1,2),
mesh(T,X,v')
xlabel('time')
ylabel('space')
zlabel('v population')
%print -dpng PopNonConst

 Nt=ceil((tf-ti)/dt);

figure,
plot(x,u(Nt/Nt,:))
hold on, plot(x,v(Nt/Nt,:))
figure,
subplot(2,2,1)
plot(x,u(0.25*Nt,:))
hold on, plot(x,v(0.25*Nt,:))
legend('population u','population v')
ylabel('population')
xlabel('space at time=3.75')


subplot(2,2,2)
plot(x,u(0.50*Nt,:))
hold on, plot(x,v(0.50*Nt,:))
legend('population u','population v')
ylabel('population')
xlabel('space at time=7.5')

subplot(2,2,3)
plot(x,u(0.75*Nt,:))
hold on, plot(x,v(0.75*Nt,:))
legend('population u','population v')
ylabel('population')
xlabel('space at time=11.25')


subplot(2,2,4)
plot(x,u(Nt,:))
hold on, plot(x,v(Nt,:))
legend('population u','population v')
ylabel('population')
xlabel('space at time=15')

%print -dpng PopNonConstAtTime

% % iteractive versione of plot
% time=1:numel(t);
% figure,
% for i=1:numel(time)
%     plot(x,u(time(i),:)), hold on
%     plot(x,v(time(i),:)),hold off
%     legend('population u','population v')
%     pause(0.0000001) %increase the value inside for more speed
%     axis([xi xf 0 2])
% end

%% plot on on space x=0.5 (useful for contant population on x)
% figure,plot(t,u(:,ceil((xf-xi)/(2*dx)))),
% hold on, plot(t,v(:,ceil((xf-xi)/(2*dx)))),
% legend('population u','population v')
% ylabel('population')
% xlabel('time')
% print -dpng PopV=0