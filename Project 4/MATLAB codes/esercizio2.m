clear all, close all
%-----------------------------------------------------------------------
Re = 1;         % Reynolds number
dt = 1e-2;      % time step
tf = 10e0;      % final time
lx = 5;         % half width of box
ly = 6;         % height of box
nx = 30;        % number of x-gridpoints
ny = 50;        % number of y-gridpoints
nsteps = 200;    % number of steps with graphic output
U0=5;
%U0=15;
%U0=40;
omega=1;
%omega=10;
kk=sqrt(omega*Re/2);
%-----------------------------------------------------------------------
nt = ceil(tf/dt); dt = tf/nt;       %aggiusto dt
x = linspace(-lx,lx,nx+1); hx = 2*lx/nx;
y = linspace(0,ly,ny+1); hy = ly/ny;
[X,Y] = meshgrid(y,x);

% soluzione esatta
uex=@(y,t) U0*exp(-kk*y).*cos(omega*t-kk*y);
%-----------------------------------------------------------------------
% initial conditions
U = zeros(nx+1,ny); V = zeros(nx,ny-1);
%dato che gli estremi est e ovest di U sono incogniti, e necessario
%aggiungerli alla nostra matrice soluzione

% boundary conditions
uN = x*0;       vN = avg(x)*0;
                vS = avg(x)*0;
uW = avg(y)*0;  vW = y*0;
uE = avg(y)*0;  vE = y*0;
%-----------------------------------------------------------------------
%condizioni al bordo per v
Vbc = dt/Re*([vS' zeros(nx,ny-3) vN']/hx^2+...
      [2*vW(2:end-1);zeros(nx-2,ny-1);2*vE(2:end-1)]/hy^2);
  
%-----------------------------------------------------------------------
%Calcolo la soluzione esatta in alcuni istanti
tt=[0 0.5 5 10];
for r= 1: length(tt)
    t=tt(r);
    ex=uex(y,t);
    
    figure
    plot(ex, y), hold on,
    axis([-lx lx 0,ly]),
    xlabel('soluzione esatta'), ylabel('y'),
    title(sprintf('t= %0.1g',t)), hold off;
  
end

pause
clear t;
%--------------------------------------------------------------------
fprintf('initialization')

%Matrix for the pressure
Lp = kron(speye(ny),K1(nx,hx,1))+kron(K1(ny,hy,1),speye(nx));

Lp(1,1) = 3/2*Lp(1,1);
perp = symamd(Lp);
Rp = chol(Lp(perp,perp));
Rpt = Rp';
Lu = speye((nx+1)*ny)+dt/Re*(kron(speye(ny),K1(nx+1,hx,1))+...
     kron(K1(ny,hy,3),speye(nx+1)));
peru = symamd(Lu);
Ru = chol(Lu(peru,peru));
Rut = Ru';

Lv = speye(nx*(ny-1))+dt/Re*(kron(speye(ny-1),K1(nx,hx,3))+...
     kron(K1(ny-1,hy,2),speye(nx)));
perv = symamd(Lv);
Rv = chol(Lv(perv,perv));
Rvt = Rv';

%vettore dell'errore
err=zeros(1,nt);
[A,B]=meshgrid(x,avg(y));

fprintf(', time loop\n--20%%--40%%--60%%--80%%-100%%\n')
for k = 1:nt
   time=(k-1)*dt;
   %la condizione al bordo varia ad ogni passo temporale
   uS = x*0+U0*cos(omega*time);
   
   % treat nonlinear terms
   Ubc = dt/Re*[2*uS' zeros(nx+1,ny-2) 2*uN']/hy^2;
   
   gamma = min(1.2*dt*max(max(max(abs(U)))/hx,max(max(abs(V)))/hy),1);
   %if gamma = 0 the we don't use UpWind 
   %if gamma = 1 we use only UpWind
   %gamma is a CFL number!
   Ue = U;
   Ue = [2*uS'-Ue(:,1) Ue 2*uN'-Ue(:,end)];
   
   Ve = [vS' V vN'];
   Ve = [2*vW-Ve(1,:); Ve; 2*vE-Ve(end,:)];
   
   Ua = avg(Ue')'; Ud = diff(Ue')'/2;
   Va = avg(Ve);   Vd = diff(Ve)/2;
   UVx = diff(Ua.*Va-gamma*abs(Ua).*Vd)/hx;
   UVy = diff((Ua.*Va-gamma*Ud.*abs(Va))')'/hy;
   Ua = avg([Ue(2,2:end-1); Ue(:,2:end-1); Ue(end-1,2:end-1)]);
   Ud = diff([Ue(2,2:end-1); Ue(:,2:end-1); Ue(end-1,2:end-1)])/2;
   Va = avg(Ve(2:end-1,:)')';
   Vd = diff(Ve(2:end-1,:)')'/2;
   U2x = diff(Ua.^2-gamma*abs(Ua).*Ud)/hx;
   V2y = diff((Va.^2-gamma*abs(Va).*Vd)')'/hy;
   %abbiamo trattato U in questo modo perche U contiene righe (agli estremi) che non
   %interagiscono con V nel calcolo della soluzione
   
   U = U-dt*(UVy+U2x);
   V = V-dt*(UVx(:,2:end-1)+V2y);
   
   % implicit viscosity
   rhs = reshape(U+Ubc,[],1);%creates and decides itself the column size
   u(peru) = Ru\(Rut\rhs(peru));%we solve the system we described before
   % but using only one line!
   U = reshape(u,nx+1,ny);
   rhs = reshape(V+Vbc,[],1);
   v(perv) = Rv\(Rvt\rhs(perv));
   V = reshape(v,nx,ny-1);
   
   % pressure correction
   rhs = reshape(diff(U)/hx+diff([vS' V vN']')'/hy,[],1);
   p(perp) = -Rp\(Rpt\rhs(perp));
   P = reshape(p,nx,ny);
   U = U-[P(2,:)-P(1,:); diff(P);P(end,:)-P(end-1,:)]/hx;
   V = V-diff(P')'/hy;
   
   % visualization
   if floor(25*k/nt)>floor(25*(k-1)/nt), fprintf('|'), end
   if k==1||floor(nsteps*k/nt)>floor(nsteps*(k-1)/nt)
     
      clf, contourf(avg(x),avg(y),P',20,'w-'), hold on
      
      Ue = [uS' avg(U')' uN'];
      Ve = [vW;avg([vS' V vN']);vE];
     % Len = sqrt(Ue.^2+Ve.^2+eps);
     % quiver(x,y,(Ue./Len)',(Ve./Len)',.4,'k')
     quiver(x,y,(Ue)',(Ve)',4,'k') %tolgo normalizzazione lunghezze
      hold off, axis equal, axis([-lx lx 0 ly])
      p = sort(p); caxis(p([8 end-7]))
      title(sprintf('t = %0.4g U = %0.2g W = %0.2g',(k-1)*dt,U0,omega))
      xlabel('x'), ylabel('y'),
      drawnow
      
   end
   
 %analizzo l'errore
 if(time>=4)
 Uex=uex(B,time);
 err(1,k)=max(max(abs(Uex'-U)));
 end
 
end
fprintf('\n')
max(err)
%=======================================================================

