%% Es 2(implementazione di tutti i metodi M1, M2, M3 e M4)
%Scommentare il metodo da eseguire insieme ad altri dati richiesti per tale
%metodo
clear all, close all

ICs=[1/2 -3];

ti=0; tf=1;

hh=[0.1 0.05 0.025 0.0125];

yEx1=@(t) (3/8)*exp(-8*t)+(1/8);
yEx2=@(t) -3*exp(-8*t);

%servono per il metodo M2
%yConv1M2=@(t) (3/8)*exp(-(8/3)*t)+1/8; 
%e' dove converge al posto di y1 il metodo M2
%yConv2M2=@(t) -3*exp(-(8/3)*t); 
%e' dove converge al posto di y2 il metodo M2

%servono per il metodo M4
% yConv1M4=@(t) (3/8)*exp(-(4)*t)+1/8;
%e' dove converge al posto di y1 il metodo M4
% yConv2M4=@(t) -3*exp(-(4)*t);
%e' dove converge al posto di y2 il metodo M4

err(:,1)=[0; hh']; %la prima colonna di errore(err) sono gli h
x=linspace(0,1);
for l=1:numel(hh)
    clear y t
    j=2; %indice per err
    h=hh(l);
    t(1)=ti;
    t(2)=ti+h;
    y(1,:)=[ICs(1) ICs(2)];
    y(2,:)=[yEx1(t(2)) yEx2(t(2))]; %questo verra riscritto durante M4
    n=ceil((tf-ti)/h);
    %% Metodo M1
%     for i=1:n-1
%         t(i+2)=t(i+1)+h;
%         y(i+2,:)=M1(@odefun,t(i),t(i+1),t(i+2),y(i,:),y(i+1,:),h);
%         if(l==1 && ((i+2)==3||(i+2)==5 || (i+2)==7 || (i+2)==9 || i+2==11))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%         if(l==2 && ((i+2)==5||(i+2)==9 || (i+2)==13 || (i+2)==17 || i+2==21))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%         if(l==3 && ((i+2)==9||(i+2)==17 || (i+2)==25 || (i+2)==33 || i+2==41))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%         if(l==4 && ((i+2)==17||(i+2)==33 || (i+2)==49 || (i+2)==65 || i+2==81))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%     end
  
    %% Metodo M2
%     clear y3;
%     y3(1,:)=ICs;%la sol di ode y'=f*3
%     y3(2,:)=[yConv1M2(t(2)) yConv2M2(t(2))];
%     for i=1:n-1
%         t(i+2)=t(i+1)+h;
%         y(i+2,:)=M2(@odefun,t(i),t(i+1),y(i,:),y(i+1,:),h);
%         %y3(i+2,:)=M2(@odefun3,t(i),t(i+1),y(i,:),y(i+1,:),h); %per ode fun*3
%         if(l==1 && ((i+2)==3||(i+2)==5 || (i+2)==7 || (i+2)==9 || i+2==11))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%         if(l==2 && ((i+2)==5||(i+2)==9 || (i+2)==13 || (i+2)==17 || i+2==21))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%         if(l==3 && ((i+2)==9||(i+2)==17 || (i+2)==25 || (i+2)==33 || i+2==41))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%         if(l==4 && ((i+2)==17||(i+2)==33 || (i+2)==49 || (i+2)==65 || i+2==81))
%             err(1,j)=t(i+2);
%             err(l+1,j)=norm([yEx1(t(i+2)) yEx2(t(i+2))]-[y(i+2,1) y(i+2,2)]);
%             j=j+1;
%         end
%     end

    %% Metodo M3
    %t(3)=t(2)+h;
    %y(3,:)=[yEx1(t(3)) yEx2(t(3))];
%     for i=1:n-2
%         t(i+3)=t(i+2)+h;
%         y(i+3,:)=M3(@odefun,t(i),t(i+2),y(i,:),y(i+1,:),y(i+2,:),h);
%         if(l==1 && ((i+3)==5 || (i+3)==7 || (i+3)==9 || (i+3)==11))%caso t=0.2 errore = 0
%             %perche e dato da yEx
%             j=j+1;
%             err(1,j)=t(i+3);
%             err(l+1,j)=norm([yEx1(t(i+3)) yEx2(t(i+3))]-[y(i+3,1) y(i+3,2)]);
%         end
%         if(l==2 && ((i+3)==5||(i+3)==9 || (i+3)==13 || (i+3)==17 || i+3==21))
%             err(1,j)=t(i+3);
%             err(l+1,j)=norm([yEx1(t(i+3)) yEx2(t(i+3))]-[y(i+3,1) y(i+3,2)]);
%             j=j+1;
%         end
%         if(l==3 && ((i+3)==9||(i+3)==17 || (i+3)==25 || (i+3)==33 || i+3==41))
%             err(1,j)=t(i+3);
%             err(l+1,j)=norm([yEx1(t(i+3)) yEx2(t(i+3))]-[y(i+3,1) y(i+3,2)]);
%             j=j+1;
%         end
%         if(l==4 && ((i+3)==17||(i+3)==33 || (i+3)==49 || (i+3)==65 || i+3==81))
%             err(1,j)=t(i+3);
%             err(l+1,j)=norm([yEx1(t(i+3)) yEx2(t(i+3))]-[y(i+3,1) y(i+3,2)]);
%             j=j+1;
%         end
%     end
    
    %% Metodo M4
%     clear y2;
%     y2(1,:)=ICs;
%     for i=1:n
%         t(i+1)=t(i)+h;
%         y(i+1,:)=M4(@odefun,t(i),y(i,:),h);
%         %y2(i+1,:)=M4(@odefun2,t(i),y(i,:),h); %la soluzione di y'=2*f;
%         if(l==1 && ((i+1)==3||(i+1)==5 || (i+1)==7 || (i+1)==9 || i+1==11))
%             err(1,j)=t(i+1);
%             err(l+1,j)=norm([yEx1(t(i+1)) yEx2(t(i+1))]-[y(i+1,1) y(i+1,2)]);
%             j=j+1;
%         end
%         if(l==2 && ((i+1)==5||(i+1)==9 || (i+1)==13 || (i+1)==17 || i+1==21))
%             err(1,j)=t(i+1);
%             err(l+1,j)=norm([yEx1(t(i+1)) yEx2(t(i+1))]-[y(i+1,1) y(i+1,2)]);
%             j=j+1;
%         end
%         if(l==3 && ((i+1)==9||(i+1)==17 || (i+1)==25 || (i+1)==33 || i+1==41))
%             err(1,j)=t(i+1);
%             err(l+1,j)=norm([yEx1(t(i+1)) yEx2(t(i+1))]-[y(i+1,1) y(i+1,2)]);
%             j=j+1;
%         end
%         if(l==4 && ((i+1)==17||(i+1)==33 || (i+1)==49 || (i+1)==65 || i+1==81))
%             err(1,j)=t(i+1);
%             err(l+1,j)=norm([yEx1(t(i+1)) yEx2(t(i+1))]-[y(i+1,1) y(i+1,2)]);
%             j=j+1;
%         end
%     end
    %% plotting
%     subplot(1,2,1),plot(t,y(:,1)),hold on
%     plot(x,yEx1(x)),hold off
%     subplot(1,2,2), plot(t,y(:,2)), hold on
%     plot(x,yEx2(x)), hold off
%     pause
end
%% stampa errore
%la prima colonna sono i vari h mentre la prima riga sono i tempi a cui e'
%calcolato l'errore
%err