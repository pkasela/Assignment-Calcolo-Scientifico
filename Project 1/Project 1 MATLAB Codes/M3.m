%% metodo M3 restituisce u_{n+3}
function yn3=M3(odefun,tn,tn2,yn,yn1,yn2,h)

yn3=-(1/4)*yn2+(1/2)*yn1+(3/4)*yn+(h/8)*[19*odefun(tn2,yn2)+5*odefun(tn,yn)];

return