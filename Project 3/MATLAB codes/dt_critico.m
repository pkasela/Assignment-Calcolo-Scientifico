function dt=dt_critico(dx,beta,D)
%questa funzione serve a calcolare il dt massimo per manenere la stabilita'

%Prende in entrata dx e beta, calcola il Pe associato, quindi il C massimo
%e infine ricava dt

 Pe=beta*dx/(2*D);
 if(Pe<0)
     dt = 0;
 end
 if(Pe<1)
     dt=Pe*dx/beta;
 else
     dt=dx/(Pe*beta);
 end
return