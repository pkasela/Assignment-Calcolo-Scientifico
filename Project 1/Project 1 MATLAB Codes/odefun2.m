%% e' la ODE che vogliamo studiare moltiplicata per 2
function f=odefun2(t,y)

f=zeros(1,2);
f(1)=2*y(2);
f(2)=2*y(2)*(y(2)-1)/(y(1));

return