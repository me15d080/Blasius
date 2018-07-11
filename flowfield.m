function[u,v,Theta]=flowfield(xo,yo,uinf,nu,correctf3,correctf5)
eta_max=yo*sqrt(uinf/xo*nu); %% similarity variable
eta_span=linspace(0,eta_max,200);
y0=[0;...
         0;...
         correctf3;...
         0;...
         correctf5];
[eta f]=rk2ode_imp2(eta_span,y0);

%% velocity components
u= f(end,2)*uinf;            
v= 0.5*(sqrt(nu*uinf/xo))*((eta_max*f(end,2))-f(end,1));
Theta=f(end,4);
end
