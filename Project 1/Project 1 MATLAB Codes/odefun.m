%% e' la ODE che studiamo
function f=odefun(t,y)

f=zeros(1,2);
f(1)=y(2);
f(2)=y(2)*(y(2)-1)/(y(1));

return