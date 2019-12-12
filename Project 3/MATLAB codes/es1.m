clear all, close all
D=1; 
%beta=0.05;
%beta=1;
beta=100;

%la soluzione esatta e':
u=@(x,t) exp(-4*pi^2*D*t).*sin(2*pi*(x-beta*t));

%determino tempo finale e quanti passi avere nel segmento (0,1)
tf=0.15;
dex=[100];%numero dei nodi

%ora utilizzo la funzione Pe_crit.m per mostrare l'area di convergenza al
%variare di Pe, C

% aux=0:0.01:10;
% for i= 1: numel(aux)
%     Pe_cr(i)=Pe_crit(aux(i));
% end
% area(aux,Pe_cr),
% xlabel('Pe'), ylabel('C'),title('area di stabilita''),
% hold off;


figure
for l=1:numel(dex)
    dx=1/dex(l);
    x=0:dx:1;
    dtcrit=dt_critico(dx,beta,D);
    
    %determino i dt da analizzare in base al valore critico
    deltat=[1.01 1 0.99 0.7 0.5 0.2 0.1 0.01]*dtcrit;
    %deltat=[1]*dtcrit;
    for k=1:numel(deltat)
        dt=deltat(k);
        t=0:dt:tf;
        Nt=numel(t);
        Nx=numel(x);
        v=zeros(Nx,Nt);
        %metto le condizioni al bordo
        v(:,1)=u(x,0)';
        v(1,:)=u(0,t)';
        v(end,:)=u(1,t)';

        for j=2:Nt
            for i=2:Nx-1
                v(i,j)=v(i-1,j-1)*(dt/(dx^2)+beta*dt/(2*dx))+v(i,j-1)*(1-2*dt/(dx^2))+v(i+1,j-1)*(dt/(dx^2)-beta*dt/(2*dx));
            end
        end

      %questa porzione di codice serve a confrontare graficamente soluzione numerica e
      %soluzione esatta 
%          figure
%          subplot(1,2,1)
         [X,T]=meshgrid(x,t);
        mesh(x,t,v'), hold on
          xlabel('x'), ylabel('t');
         %subplot(1,2,2)
%          mesh(X,T,u(X,T))
%          figure
%          mesh(x,t,v'), hold on
%         mesh(X,T,u(X,T))


        %calcolo ora %l'errore tra sol esatta e approssimazione
        A=u(X,T)';
        B=A-v;
        err(l,k)=max(max(abs(B)));
    end
    
    %costruisco ora il grafico che confronta i vari dt con l'errore
    str_title=strcat('delta x=',sprintf('%1.0g',dx));
    subplot(2,ceil(numel(dex)/2),l)
    loglog(deltat,err(l,:),'*-')
     title(str_title)
     xlabel('delta t'),ylabel('err')
    
end
