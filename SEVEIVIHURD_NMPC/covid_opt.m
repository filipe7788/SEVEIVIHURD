function cost=covid_opt(x,t,Model_0,N0,Usig,par,N_obd,Q)
global psir costd

for k=1:size(Usig,2)

    x=Usig(k);

[t,Psi,S,E,V,Ev,Ia,Is,Iav,Isv,H,U,R,Rv,D,Nw,NwV]=solver_model(x,t,Model_0,N0,par,N_obd);

if t(1) < 70
nlhos = 466;
nluti = 422;

for i=1:1:length(H)

if (H(i)>nlhos)||(U(i)>nluti)
    delta1(i) = Is'*Is + Isv'*Isv + H'*H + U'*U;  % Include Vaccinated Infected
else
    delta1(i)=0;
end
end

cost = Is(end) + Isv(end) + Q(1)*x + sum(delta1);

elseif t(1)<170  %214   
    
nlhos = 1610;
nluti = 1210;

for i=1:1:length(H)

if (H(i)>nlhos)||(U(i)>nluti)
     delta2(i) = Is'*Is + Isv'*Isv + H'*H + U'*U;  % Include Vaccinated Infected
else
    delta2(i)=0;
end
end

        cost = Is(end) +  Isv(end) + Q(2)*x + sum(delta2);

else

    nlhos = 1610;
    nluti = 1210;
    
    for i=1:1:length(H)

        if (H(i)>nlhos)||(U(i)>nluti)
             delta2(i) = Is'*Is + Isv'*Isv + H'*H + U'*U;  % Include Vaccinated Infected
        else
            delta2(i)=0;
        end
        end

       cost = Is(end) + Isv(end) + Q(3)*x + sum(delta2);

end

if cost < costd
    psir = x(1);
    costd=cost;
end


end
end