%% Solver Model SIRD
function [t,Psi,S,E,V,Ev,Ia,Is,Iav,Isv,H,U,R,Rv,D,Nw,NwV]=solver_model(x,t,Model_0,N0,par,N_obd)


% Resolve o modelo de no tempo t
%options=odeset('NonNegative',(1:8));
[t,y] = ode45(@(t,y)covid_odes(t,y,x,N0,par,N_obd),t,Model_0);

Psi=y(:,1);
S=y(:,2);
E=y(:,3);
V=y(:,4);
Ev=y(:,5);
Ia=y(:,6);
Is=y(:,7);
Iav=y(:,8);
Isv=y(:,9);
H=y(:,10);
U=y(:,11);
R=y(:,12);
Rv=y(:,13);
D=y(:,14);
Nw=y(:,15);
NwV=y(:,16);

end
