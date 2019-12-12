function C = Pe_crit(Pe)

%questa funzione calcola il C massimo per mantenere la stabilita' dato
%il Pe in entrata

if(Pe<=0)
    C=0;
end
if(Pe<=1)
    C=Pe;
else
    C=1./Pe;
end

