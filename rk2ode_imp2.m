function[eta,f]=rk2ode_imp2(eta_span,y0)
%%  Integrator based on Heun's method for eq 1.2 (in report)
global Pr
f1=0*eta_span;
f2=f1;
f3=f1; 
f4=f1;
f5=f1;

%% Initial Value
f1(1)=y0(1);
f2(1)=y0(2);
f3(1)=y0(3);
f4(1)=y0(4);
f5(1)=y0(5);

N=length(eta_span);
deta=eta_span(2)-eta_span(1);

for i=1:N-1
%% Predictor step  
%% Guess the value of function at next point     
gf1=f1(i)+f2(i)*deta;
gf2=f2(i)+f3(i)*deta;
gf3=f3(i)*(1-f1(i)*0.5*deta);
gf4=f4(i)+f5(i)*deta;
gf5=f5(i)*(1-0.5*Pr*f1(i)*f5(i))*deta;

%% Slope at Guessed point 
gf1prime=gf2; 
gf2prime=gf3;
gf3prime=-0.5*gf1*gf3;
gf4prime=gf5;
gf5prime=-0.5*Pr*gf1*gf5;

%% Corrector step
%% Average the slope at two ends of interval 

    cf1prime=0.5*(gf1prime+f2(i));
    cf2prime=0.5*(gf2prime+f3(i));
    cf3prime=0.5*(gf3prime-0.5*f1(i)*f3(i));
    cf4prime=0.5*(gf4prime+f5(i));
    cf5prime=0.5*(gf5prime-0.5*Pr*f1(i)*f5(i));

%% Calculate the value at next step
 f1(i+1)= f1(i)+cf1prime*deta;
 f2(i+1)= f2(i)+cf2prime*deta;
 f3(i+1)= f3(i)+cf3prime*deta;
 f4(i+1)= f4(i)+cf4prime*deta;
 f5(i+1)= f5(i)+cf5prime*deta;  

 
end
eta=eta_span;
f=[f1' f2' f3' f4' f5'];
end
