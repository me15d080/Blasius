function correctf3=get_correctf3( guess1,guess2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% module for using eq 1.1 (in report)for getting %        
% correct f" value by shooting method            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eta_span=linspace(0,10,500);          % Domain of integration
fprime0=[0;...                        % First Initial guess 
    0;...
    guess1];

[eta,f]=rk2ode_imp(eta_span,fprime0); % Integrate over domain of integration
yl=1-f(end,2);                        % See whether we are ahead/ behind the correct 

fprime0=[0;...                        % Second Guess
    0;...
    guess2];                   
[eta,f]=rk2ode_imp(eta_span,fprime0); % Integrate over domain of integration
yu=1-f(end,2);                        % See whether we are ahead/ behind the correct

%% Check whether root is bracketed or not
if (yl*yu>0)
    fprintf('root not bracketed\n')
    
    return
end

%% Root has been bracketed:Go for Bisection method
iter=0;
err=abs(guess1-guess2);
while(err>1e-5)
    corrf3=0.5*(guess1+guess2); %% Next guess "bisects" the interval  
    fprime0=[0;...              %% Next IC
        0;...
        corrf3];
    
    [eta,f]=rk2ode_imp(eta_span,fprime0);% Integrate over domain of integration
    yr=1-f(end,2);
%% Which part contains the root ? (guess1, corrf3) or (corrf3, guess2) 
    if (yl*yr>=0)
        guess1=corrf3;
        yl=yr;
    else
        guess2=corrf3;
        yu=yr;
    end
%% Check whether we need to do next iteration
    err=abs(guess1-guess2); 
    iter=iter+1;
end
%% Solution
correctf3=corrf3;
end
