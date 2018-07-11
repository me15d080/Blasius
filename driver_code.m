%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Driver code:                                                                         %
% This is the main solver which calculates flow and temp field. Please refer           %
% to the report for symbols, governing equations and numerical strategy.               %
%                                                                                      %
%  Module/Variable         Description                                                 %
%  Modules :                                                                           %

%  rk2ode_imp.m            Integrator based on Heun's method for eq 1.1 (in report)    %
%  rk2ode_imp2.m           Integrator based on Heun's method for eq 1.2 (in report)    %
%  get_correctf3.m         module for using eq 1.1 (in report)for getting              %
%                          correct f" value by shooting method                         %
%  get_correctf5.m         module for using eq 1.2 (in report)for getting              %
%                          correct f" value by shooting method                         %
%  flowfield.m             calculates the velocity components and temperature at each  %
%                          grid point in computational domain                          %
%  usr_input.m             Contains all the inputs needed from user                    %

%   Variables:
%  guess1_f3, guess2_f3    refer to usr_input                                          % 
%  guess1_f5, guess2_f5    refer to usr_input                                          %
%  uinf,Tinf,nu,Lx,Ly,Nx,Nyrefer to usr_input                                          %
%  f                       refer to report eq 1.1 and 1.2                              %
%  correctf3               correct f"(0) for eq 1.1 (in report)                        %
%  correctf5               correct theta'(0) for eq 1.2 (in report)                    %
%  eta                     similarity variable                                         %
%  eta_span                Domain of integration                                       %
%  yo                      Initial conditions                                          %
%  Cf_x                    Local skin friction coefficient                             %
%  Nu_x                    Local Nusselt number                                        %
%  sf                      stream function                                             %
%  U,V                     x-y Velocity components                                     %
%  T ,Tplate               Temp @ grid point and @ surface respectively                %
%  Theta                   Non dimensional temperature                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all;close all;

%% Read the user input

usr_input 

%% Shooting for correct f(3) 

correctf3=get_correctf3(guess1_f3,guess2_f3);

%% Shooting correct f(5)

correctf5=get_correctf5 (guess1_f5,guess2_f5,correctf3);

%% Calculate non dimensional velocity and temp profile
%% 
yo   =  [0;...
         0;...
         correctf3;...
         0;...
         correctf5];                                         %% correct initial condition

     
eta_span=linspace(0,15,500);                                 %% Domain of integrtaion

[eta, f]=rk2ode_imp2(eta_span,yo);                           %% Integrate to solve eq 1.2

%% Post Processing 1:Non dimensional velocity and temp profile on the same plot
figure(1)
plot(f(:,2),eta);                                            %% f(:,2) stores f2 i.e. Non dimensional velocity u/uinf
axis([0 1.2 0 15])
hold on
plot(f(:,4),eta,'g');                                        %% f(:,2) stores f2 i.e. Non dimensional velocity u/uinf
xlabel('Non dimensional Temp');ylabel('eta (similarity variable)');
title('Eta Vs Non dimensional velocity')
hold off

%% Create a 2d computational grid in physical space

x=linspace(0.5,Lx,Nx);
y=linspace(Ly,0,Ny);
[X,Y]=meshgrid(x,y);

%% Initialise variables (Refer Description at the top for meaning of variables)
U=X;
V=U;
T=U;
Cf_x=zeros(1,Nx);
Nu_x=Cf_x;
h_x=Nu_x;
sf=U;
X_s=Cf_x;

%% Calculate Flow field, Local skin friction coefficient and Local Nusselt Number

for i=1:Nx

   X_s(i)=X(end,i);
   Re_x=(uinf*X_s(i))/nu;
   Cf_x(i)=2*correctf3/sqrt(Re_x);
   Nu_x(i)=correctf5*sqrt(Re_x);
   h_x(i)=k*Nu_x(i)/X_s(i);

         for j=1:Ny

             x0=X(i,j);y0=Y(i,j);
             [u,v,Theta]=flowfield(x0,y0,uinf,nu,correctf3,correctf5);
             U(i,j)=u;
             V(i,j)=v;
             T(i,j)=Tplate+(Tinf-Tplate)*Theta;
             
         end

end

%% Post Processing 2: Cf_x and Nu_x
figure(2)
plot(X_s,Cf_x);
xlabel('x');ylabel('Cf_x');
title('Variation of local skin friction coefficient with x');
 
figure(3)
plot(X_s,Nu_x);
xlabel('x');ylabel('Nu_x');
title('Variation of local Nusselt Number with x');

figure(4)
plot(X_s,h_x);
xlabel('x');ylabel('h_x');
title('Variation of heat transfer coefficient with x');
 


%% Post Processing 3:Velocity and Temperature contour plots
 figure(5)
 contourf(X,Y,U,25)                 %% Velocity contours for determining boundary layer thickness
 xlabel('X');ylabel('Y');title('Velocity contours')
 figure(6) 
 contourf(X,Y,T,25)                 %% Temperature contours
 xlabel('X');ylabel('Y');title('Temperature contours')

%% For streamlines (BONUS PROBLEM)
 yi=[0;...
     0;...
     correctf3];
 for i=1:Nx
 
          for j=1:Ny
 
               x0=X(i,j);y0=Y(i,j);
               
               eta0=y0*sqrt(uinf/(nu*x0));
               eta_span=linspace(0,eta0,200);
               [eta, f]=rk2ode_imp(eta_span,yi);
               sf(i,j)=f(end,1)*sqrt(nu*uinf*x0);
              
          end
 
 end
%% Post Processing 4: streamlines

 figure (7)
 contourf(X,Y,sf,25)%% streamlines
 xlabel('X');ylabel('Y');title('streamlines')
 
%% Save the Data
% save('U_Tdata','U','T')



