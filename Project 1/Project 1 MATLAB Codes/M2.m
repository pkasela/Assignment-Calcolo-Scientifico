%% metodo M2 restituisce u_{n+2}
function yn2=M2(odefun,tn,tn1,yn,yn1,h)

yn2=yn1+(h/3)*(3*odefun(tn1,yn1)-2*odefun(tn,yn));

return