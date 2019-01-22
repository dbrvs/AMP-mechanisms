% function that calculates the SILVE odes 
function dX=infection_model(~,X,...
    aS,dS,bt,rho,tau,dI,k,xi,aL,dL,w,aE,dE,E50,p,g)

    S   = X(1); %hiv specific cd4+/ccr5+ T cells
    I_U = X(2); %infected unproductive viral strain
    I_P = X(3); %infected productive strain
    L_U = X(4); %latently infected unproductive strain
    L_P = X(5); %latently infected productive strain
    E   = X(6); %adaptive immune response
    V   = X(7); %virus
  
dX=zeros(7,1);

dX(1) = aS - dS*S - bt*S*V;% + xi*E; %dS
dX(2) = (1-rho)*tau*bt*S*V - dI*I_U - k*I_U*E + xi*L_U*E; %dI_U
dX(3) = (1-rho)*(1-tau)*bt*S*V - dI*I_P - k*I_P*E + xi*L_P*E; %dI_P
dX(4) = rho*tau*bt*S*V + aL*L_U - dL*L_U - xi*L_U*E; %dL_U
dX(5) = rho*(1-tau)*bt*S*V + aL*L_P - dL*L_P - xi*L_P*E; %dL_U
dX(6) = aE - dE*E + w*E*(I_U+I_P)/(E+E50);
dX(7) = p*I_P - g*V - bt*S*V;
       
    
