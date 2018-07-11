% file name : usr_input.m                                                       
% Description of code:                                                          
% This file is meant to specify/ change the inputs to be provided               
% by the user. All changes are to be made here only and no where else           
% in the rest of code. Once user is satisfied with the changes,he should        
% run the driver_code.m for the numerical solution of boundary layer equations.   

% Geometry
Lx=10;              %length of computational domain
Ly=4;               % height of computational domain
Nx=100;             % No. of x grid points
Ny=100;             % No. of y grid points
% Flow property
uinf=35;            % free stream  velocity
Tinf=50;            % free stream temp
Tplate=30;          % flat plate surface temperature
% Fluid property
nu=0.8 ;            % kinematic viscosity 
k=1;                % Fluid thermal conductivity
global Pr;          % Prandtl Number (k/rho*cp)
Pr=0.1;
% Guesses for Shooting Method %%

% User guess for f double prime

guess1_f3=0;
guess2_f3=5;

% User guess for Theta prime

guess1_f5=0;
guess2_f5=5;

% Note 
% f means non dimensional stream function used in Blasius momentum equation.
% theta means non dimensional temperature used in Pohlhausen Energy equation.
% f3 mean (d2f/dn2) and f5 means (d2(theta)/d(eta2))
