%solve the infection model
%DBR 6/2017


%model parameters
aS=100;   %S birth rate
dS=0.3;  %S death rate
bt=1e-4;  %infectivity 
rho=1e-4;  %fraction latent
tau=0.05; %fraction productive
dI=1;
k=1;
xi=1e-5;
aL=0.015;
dL=0.0155;
w=1;
aE=1e-5;
dE=0.003;
E50=250;
p=3e3;
g=23;

ops = odeset('RelTol', 1e-10); %tolerances for ODE solver

tspan=linspace(0,30,1e3);

I0=1e-3; %1 infected cell per mL
X0=[aS/dS,0,I0,0,0,aE/dE,0];
    
[~, X]=ode23s(@(t,X)infection_model(...
        t,X,aS,dS,bt,rho,tau,dI,k,xi,aL,dL,w,aE,dE,E50,p,g),tspan,X0,ops);

    
subplot(231)
semilogy(tspan,X(:,1)) %S
title('CD4+/CCR5+ HIV specific T cells')

subplot(232)
semilogy(tspan,X(:,2),tspan,X(:,3)) %S
title('infected cells')
legend('unproductive','productive')

subplot(233)
semilogy(tspan,X(:,4),tspan,X(:,5)) %S
title('latently infected cells')
legend('unproductive','productive')

subplot(234)
semilogy(tspan,X(:,6)) %S
title('adaptive immune cells')

subplot(235)
semilogy(tspan,X(:,7)*1e3) %S
title('viral RNA')



