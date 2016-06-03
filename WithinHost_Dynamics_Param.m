function F = WithinHost_Dynamics_Param(k_p, a_p, rL_p, start, duration, ART, initialStrain, drugM)

% Solves within-host dynamics for 2 strains resistant and not resistant to
% drug

%%% PARAMETERS %%%

ns = 2;             % number of strains
startstrain = initialStrain;    % Strain from which the infection is started
ntsteps = 1e4;         % number of time steps

yfact = 365;           % factor to convert from per day to per year

forward_mut = (5e-5)^1;    % Mutation rate of mutating to the next strain (index 1 higher)
back_mut = (5e-5)^1;       % Mutation rate of mutating into a strain with index 1 lower (NOTE: For non-neutral evolution we set forward_mut = back_mut)
rand_mut = 0;          % Small mutation rate, mutation rate of mutating directly to a random other strain

drug = @(t) drugM*concentration(t,start,duration,48/8760,0,1.0); % The concentration of PrEP
g = @(t) yfact * [1.0*(1-drug(t)); 0.9]; % The fitness of strain 2 does not depend on PrEP (it is resistant)

concentration(t,Delay,inf_test,48/8760,0,1.0)
%%% Reservoir parameters %%%

% Probability that newly infected cells become latent, vector for three stages (acute, asymptomatic, AIDS). [Note: in this work we do not change the values of parameters over the three stages]                  
k = k_p*[1,1,1];
% Activation rate of latent cells (per day), vector for three stages (acute, asymptomatic, AIDS).
a = a_p*[1,1,1];
% Size ratio between latent and active compartment: (size of latent reservoir)/(size of active compartment)
rL = rL_p*[1,1,1];


%%% Duration of infection (as determined by spVL) %%%

D1 = 0.25;                                          % Duration of acute infection (years)
D2 = 4;%Hill_func_down(10.^spvl,25.4,3058,0.41);       % Vector of durations of asymptomatic phase (years)
D3 = 0.75;                                          % Duration of AIDS (years)


%%% Time vectors %%%

T1 = D1*ones(ns,1);             % Time of acute infection (vector)
T2 = T1 + D2;                   % Time of acute inf + asymp phase (vector)
T3 = T2 + D3*ones(ns,1);        % Time of complete infection (vector)

Tmax = ART;             % Maximum time, after this time we stop all calculations
dtau = Tmax/ntsteps;            % time step
tau = 0:dtau:Tmax;              % time vector, with all time steps
ltau = length(tau);             % length of time vector

tau1 = round(ntsteps*T1/Tmax);        % Index of time acute -> asymptomatic
tau2 = round(ntsteps*T2/Tmax);        % Index of time asymptomatic -> AIDS
tau3 = round(ntsteps*T3/Tmax);        % Index of time AIDS -> death




%%% MODEL %%%


% Define inital conditions: start with single virion of certain strain in
% both compartments (active comp = 1:ns, reservoir = (ns+1):2*ns).
Init_temp = eye(ns);
Init = zeros(2*ns,ns);
Init(1:ns, :) = Init_temp;
Init((ns+1):2*ns, :) = Init_temp;

% Solve the within-host ODEs with ODE solver
x = zeros(ns, ltau, ns);
y = zeros(ns, ltau, ns);
tempstart = zeros(2*ns,1);

   Times = 1:min(tau1(startstrain),length(tau)); % We want the solution to run until ART
   [x(:,Times,startstrain), y(:,Times,startstrain)] = solve_qsODE_PrEP(tau(Times),Init(:,startstrain),ns,g,forward_mut,back_mut,rand_mut,k(1),yfact*a(1),rL(1));       % Acute phase
   
   if tau1(startstrain)+1 <= length(tau)
   tempstart(1:ns) = x(:,tau1(startstrain),startstrain);
   tempstart((ns+1):2*ns) = y(:,tau1(startstrain),startstrain);
   Times = (tau1(startstrain)+1):min(tau2(startstrain),length(tau)); 
   [x(:,Times,startstrain), y(:,Times,startstrain)] = solve_qsODE_PrEP(tau(Times),tempstart,ns,g,forward_mut,back_mut,rand_mut,k(2),yfact*a(2),rL(2));       % Acute phase
   end
   
   if tau2(startstrain)+1 <= length(tau)
   tempstart(1:ns) = x(:,tau2(startstrain),startstrain);
   tempstart((ns+1):2*ns) = y(:,tau2(startstrain),startstrain);
   Times = tau2(startstrain)+1:min(tau3(startstrain),length(tau));
   [x(:,Times,startstrain), y(:,Times,startstrain)] = solve_qsODE_PrEP(tau(Times),tempstart,ns,g,forward_mut,back_mut,rand_mut,k(3),yfact*a(3),rL(3));       % Acute phase
   end

   
% figure
% plot(tau,drug(tau));
% 
% figure
% plot(tau,x(:,1:end,startstrain));
% 
% 
% figure
% plot(tau,y(:,1:end,startstrain));

F = [x(2,end,startstrain) y(2,end,startstrain)]; % Output the frequency of strain 2 at time ART when initially infected with strain 1

end