%% metodo M1 restituisce u_{n+2}
function yn2=M1(odefun,tn,tn1,tn2,yn,yn1,h)

%usa fsolve per cercare u_{n+2}
options = optimset('Display','off'); %opzione per disabiltiare le scritte dovute a fsolve
fun=@(x) x+yn1-2*yn-(h/4)*(odefun(tn2,x)+8*odefun(tn1,yn1)+3*odefun(tn,yn)); %la funzione di cui
%calcolare lo zero per noi x e' u_{n+2}
yn2=fsolve(fun,yn,options);%uso questo metodo perch? il metodo di Newton non funzionava
%con l'errore della matrice Jacobiana singolare per h=0.025 e 0.0125

return