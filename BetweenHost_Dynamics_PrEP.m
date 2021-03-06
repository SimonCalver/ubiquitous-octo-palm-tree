%%% FIGURE 5 %%%
% Code adapted from Lythgoe, Pellis and Fraser (2013), "Is HIV
% short-sighted? Insights from a multistrain nested model.".

% Calculate the strain distribution in a host population

% close all;
clear all;

%%% PARAMETERS %%%

ns = 2;            % number of strains
s = (1:ns)';        % Column vector containing the strain index
startstrain = 1;    % Strain from which the infection is started

gmin = 1;              % growth rate of the least fit strain (per day)
gmax = 1.05;           % growth rate of the fittest strain (per day)
yfact = 365;           % factor to convert from per day to per year
spvlmin = 2;           % minimum log(spvl), assume that these are all the same so we only worry about resisitance to drug, (i.e max = min)?
spvlmax = 7;           % maximum log(spvl)

forward_mut = 5e-5;    % Mutation rate of mutating to the next strain (index 1 higher)
back_mut = 5e-5;       % Mutation rate of mutating into a strain with index 1 lower (NOTE: For non-neutral evolution we set forward_mut = back_mut)
rand_mut = 0;          % Small mutation rate, mutation rate of mutating directly to a random other strain

% g = yfact * vir_repr_rate_index_lin(s,gmin,gmax); % Wihtin-host strain fitness, needs to be function of time (maybe)

drug = @(t) concentration(t,0.02,0.8,0,1);% 0; %concentration(t,0.2,15,0,0.05);

%1;% heaviside(2-t) - heaviside(1-t) + heaviside(4-t) - heaviside(3-t) + heaviside(6-t) - heaviside(5-t) + heaviside(8-t) - heaviside(7-t) + heaviside(10-t) - heaviside(9-t) ...
%           + heaviside(12-t) - heaviside(11-t) + heaviside(14-t) - heaviside(13-t) + heaviside(16-t) - heaviside(15-t) + heaviside(18-t) - heaviside(17-t) + heaviside(20-t) - heaviside(19-t); %min(t/10,1);%

% abs(sin(t/400));% abs(sin(t/(yfact*5))); %1 - 200./(200+t.^2); % the varying concentration of drug
g = @(t) yfact * [1.05*(1-drug(t)); 1];
% max( yfact * (vir_repr_rate_index_lin(s,gmin,gmax)*(1-drug(t)) + vir_repr_rate_index_lin(s,gmin,0.95)*drug(t)), 0 ); 
% .* ((1-(s-1)./(ns-1)) - drug(t)) ),0); 
% (1 - drug(t)),0); %(350 + ( ns - 1 ) / ( ns - 1 ) * ( s - 1 ))*concentration(t),0); % the reproduction rate is decreased by the drug, more for 
                                                                                     % some strains than others


spvl = spvl_index_lin(s,spvlmin,spvlmax); % Strain set-point viral loads

N0 = 10000; % Total initial population size
B = 200; % Constant total birth rate per year
nu = 0.02; % Constant per capita yearly death rate

dt = 0.01; % Time-step for the computation of the full dynamics
tmin = 0;
tmax = 300;


%%% Reservoir parameters %%%

% Probability that newly infected cells become latent, vector for three stages (acute, asymptomatic, AIDS) [Note: in this work we do not change the values of parameters over the three stages]                  
k = 1e-4*[1,1,1];
% Activation rate of latent cells (per day), vector for three stages (acute, asymptomatic, AIDS)
a = 1e-2*[1,1,1];
% Size ratio between latent and active compartment: (size of latent compartment)/(size of active compartment)
rL = 0.1*[1,1,1];


%%% Predefined infectivity profile (alphas),  assume that these are all the same %%%

D1 = 0.25;                                          % Duration of acute infection (years)
D2 = Hill_func_down(10.^spvl,25.4,3058,0.41);       % Vector of durations of asymptomatic phase (years)
D3 = 0.75;                                          % Duration of AIDS (years)

lambda1 = 2.76;                                     % Infectiousness during acute stage
lambda2 = Hill_func_up(10.^spvl,0.317,13938,1.02);  % Vector of infectiousness during asymptomatic phase
lambda3 = 0.76;                                     % Infectiousness during AIDS

% tot = tic;

%%% WITHIN-HOST DYNAMICS %%%

%%% Time vectors %%%

T1 = D1*ones(ns,1);
T2 = D1 + D2;
T3 = D1 + D2 + D3*ones(ns,1);

appT1 = round(T1/dt) * dt; % Approximate time of T1 using the coarser grid 
appT2 = round(T2/dt) * dt;
appT3 = round(T3/dt) * dt;
appTmax = max(appT3);
tt = 0:dt:appTmax;
ltt = length(tt);
tau1 = round( appT1 / dt )+1; % index of T1 in the grid
tau2 = round( appT2 / dt )+1; % index of AT2 in the grid
tau3 = round( appT3 / dt )+1; % index of AT3 in the grid
% I have put the "+1"s in because there should be 26 timesteps (when dt =
% 0.01) for D1 but this only gives 25
disc_vec = exp(-nu*tt);    % natural mortality over the course of an infection
disc_last = exp(-nu*appT3); % The last point of the vector of discount factor due to natural mortality
alphas = zeros(ns,ltt); % Matrix, with overall time varying infectivity profiles, one line per each SPVL
for i = 1:ns
    alphas(i,1:tau1(i)) = lambda1;      % I chose that the step function is right continuous (the value in an interval is that at the right end). The acute phase has no area, if lT1=1
    alphas(i,(tau1(i)+1):tau2(i)) = lambda2(i);
    alphas(i,(tau2(i)+1):tau3(i)) = lambda3;
    alphas(i,:) = alphas(i,:) .* disc_vec;      % Incorporate natural mortality in total infectivity, because these terms always occur together
end

% tic;

% Initial conditions - start with 1 infecting virion (in each compartment)
Init_temp = eye(ns);
Init = zeros(2*ns,ns);
Init(1:ns, :) = Init_temp;
Init((ns+1):2*ns, :) = Init_temp;

% Solve within-host system, using ODE solver

x = zeros(ns,ltt,ns); % One 2D matrix for each starting point
y = zeros(ns,ltt,ns); % One 2D matrix for each starting point
tempstart = zeros(2*ns,1);

for i = 1:ns
   [x(:,1:tau1(i),i), y(:,1:tau1(i),i)] = solve_qsODE_PrEP(tt(1:tau1(i)),Init(:,i),ns,g,forward_mut,back_mut,rand_mut,k(1),yfact*a(1),rL(1));       % Acute phase
   tempstart(1:ns) = x(:,tau1(i),i);
   tempstart((ns+1):2*ns) = y(:,tau1(i),i);
   [x(:,(tau1(i)+1):tau2(i),i), y(:,(tau1(i)+1):tau2(i),i)] = solve_qsODE_PrEP(tt((tau1(i)+1):tau2(i)),tempstart,ns,g,forward_mut,back_mut,rand_mut,k(2),yfact*a(2),rL(2));       % Chronic phase
   tempstart(1:ns) = x(:,tau2(i),i);
   tempstart((ns+1):2*ns) = y(:,tau2(i),i);
   [x(:,tau2(i)+1:tau3(i),i), y(:,tau2(i)+1:tau3(i),i)] = solve_qsODE_PrEP(tt((tau2(i)+1):tau3(i)),tempstart,ns,g,forward_mut,back_mut,rand_mut,k(3),yfact*a(3),rL(3));       % AIDS phase
end

% time_ODE = toc;
% disp(['Time to solve within-host system:          ',num2str(time_ODE,'%.2f'),' seconds.'])
% disp(' ');

% Determine strain-specific infectivity profiles and next-generation matrix
betas = zeros(ns,ltt,ns);   % strain-specific infectivity profiles (with natural mortality incorporated)
Kc = zeros(ns);             % Next generation matrix
for i = 1:ns
    for j = 1:ns
        for itt = 1:tau3(j)
            betas(i,itt,j) = alphas(j,itt) * x(i,itt,j) ; 
        end
        Kc(i,j) = dt * trapz( betas(i,1:tau3(j),j) );
    end
end

%%% BETWEEN-HOST DYNAMICS %%%

%%% Initialisation %%%

ICint = appTmax; % Duration of interval of initial conditions, before tmin, during which incidence is constant at I0
iT0 = ICint/dt;

t = (tmin-ICint):dt:tmax;
lt = length(t);
I = zeros(ns,lt);       % Matrix of incidence over absolute time, one row per type
Ibar = zeros(ns,lt);    % Matrix of relative incidence over time
Itot = zeros(1,lt);     % Total incidence
J = zeros(ns,lt);       % Prevalence is 0 at the beginning, because there are only newly infected initial infective
Jbar = zeros(ns,lt);    % Relative prevalence
Jtot = zeros(1,lt);     % Total prevalence
S = zeros(1,lt);        % Total susceptibles are N at the beginning, for the same reason
Y = zeros(ns,lt);       % Total amount of virus in the population, one row per strain
N = zeros(1,lt);        % Total population size, which is variable and starts at N0

%%%initial conditions %%%

% tic;

% Basic initial conditions
I0 = zeros(ns,1);   % Column vector of initial values for each strain
I0(startstrain) = 1;          % Initial condition for the incidence in the population, with only strain 1 present
S0 = N0;
for iic = 1:iT0
    I(:,iic) = I0;
    Ibar(:,iic) = I0 / sum(I0);
    Itot(iic) = sum(I0);
    N(iic) = N0;
    S(iic) = S0;
    J(:,iic) = zeros(ns,1);         % Could do better, probably
    Jbar(:,iic) = zeros(ns,1); 
    Jtot(iic) = 0; 
end
N(iT0) = N0;
S(iT0) = N(iT0) - Itot(iT0) * dt * trapz( disc_vec );

% Solve between-host system
FOI = zeros(ns,lt);         % Force of infection of different strains over time
FOI_temp2 = zeros(ns,1);
for it = (iT0+1):lt
    discI = zeros(ns,1);
    for i = 1:ns
        FOI_temp = zeros(ns,ltt);       % FOI from type j at any point in the past (acting on i, but not depending explicitly, i.e. overwritten every loop)
        for j = 1:ns
            FOI_temp(j,1) = 0;          % The first point of this vector is always 0, because it refers to the time point t(it), the one we are calculating the solution at
            for itt = 2:min(it,tau3(j))
                FOI_temp(j,itt) = betas(i,itt,j) * I(j,it-itt+1);
            end
            FOI_temp2(j) = dt * trapz( FOI_temp(j,1:min(it,tau3(j))) );     % Contibution of j to the FOI acting on i
        end
        FOI(i,it) = sum( FOI_temp2 );
        I(i,it) = S(it-1) * FOI(i,it) / N(it-1);
        J(i,it) = dt * trapz( I(i,max(1,it:-1:it-tau3(i)+1)) .* disc_vec(1:tau3(i)) );
        discI(i) = I(i,it-tau3(i)+1) * disc_last(i);% beause of the "+1"s I added earlier I need to take 1 from tau3 to make it the same. 
        % This may well be wrong. Thse lines used to look like this:
        % J(i,it) = dt * trapz( I(i,max(1,it:-1:it-tau3(i))) .* disc_vec(1:tau3(i)+1) ); % I suspect the max is not necessary
        % discI(i) = I(i,it-tau3(i)) * disc_last(i);
    end
    Itot(it) = sum(I(:,it));
    Jtot(it) = sum(J(:,it));
    S(it) = N(it-1) - Jtot(it);
    N(it) = N(it-1) + dt * ( B - nu * N(it-1) - sum( discI ) );
    Ibar(:,it) = I(:,it) / Itot(it);
    Jbar(:,it) = J(:,it) / Jtot(it);
%     disp(['Step ',num2str(it-iT0),' of ',num2str(lt-iT0),' completed'])
end
% time_dyn = toc;

% Calculate R0 and the average spVL at the end of the simulation
R0_end = N(lt)/S(lt);
avSPVL_end = (J(:,lt)'*spvl)/sum(J(:,lt));

% disp(' ');
% disp(['Time to solve the full dynamics:             ',num2str(time_dyn,'%.2f'),' seconds.'])
% disp(' ');

%%% FIGURES %%%

%{
figure
clf;
hold on;
plot(t,N/1000,'k','Linewidth',2);
plot(t,Jtot/1000,'k--','Linewidth',2);
axis([ 0 tmax 0 10.05]);
set(gca,'XTick',0:25:tmax,'YTick',0:2:10,'box','off','Fontsize',12);
title('Population size and prevalence','Fontsize',20);
xlabel('Time (years)','Fontsize',16);
ylabel('Thousands','Fontsize',16);
legend('Total population','Total prevalence','Location','NorthEast')
hold off;
%}



scrsz = get(groot,'ScreenSize');
main = figure('Position',[1 0 scrsz(3) scrsz(4)]);



cfig = subplot(2,2,1);
set(gcf,'DefaultAxesColorOrder',jet(ns))
plot(t,Jbar,'Linewidth',2.5);
axis([ 0 tmax 0 1.005]);
set(gca,'XTick',0:25:tmax,'YTick',0:0.2:1,'box','off','Fontsize',12);
title('Stratified prevalence','Fontsize',20);
xlabel('Time (years)','Fontsize',16);
ylabel('Strain frequency','Fontsize',16);
% Legend
legendentries = cell(ns,1);
for i = 1:ns
    legendentries(i,1) = { ['Strain ',num2str(i)] };
end
legend(legendentries,'Location','EastOutside','FontSize',12);
str(1) = {['R0 = ',num2str(R0_end)]};
str(2) = {['log(spVL) = ',num2str(avSPVL_end)]};
R0text = text(80,0.6,str,'FontSize',14);
set(R0text,'EdgeColor',[0 0 0])


fitness = zeros(ns,length(tt));

for i = 1:length(tt)
    rep = g(tt(i));
    for k = 1:length(s)
    fitness(k,i) = rep(k);
    end
end


ffig=subplot(2,2,2);
set(gcf,'DefaultAxesColorOrder',jet(ns))
plot(tt,fitness,'Linewidth',2.5);
axis([ 0 inf -inf inf]);
set(gca,'box','off','Fontsize',12);
title('','Fontsize',20);
xlabel('Time (years)','Fontsize',16);
ylabel('Fitness','Fontsize',16);
% Legend
legendentries = cell(ns,1);
for i = 1:ns
    legendentries(i,1) = { ['Strain ',num2str(i)] };
end
legend(legendentries,'Location','EastOutside','FontSize',12);



ifig = subplot(2,2,3);
set(gcf,'DefaultAxesColorOrder',jet(ns))
plot(tt,x(:,:,1),'Linewidth',2.5);
axis([ 0 inf -inf inf]);
set(gca,'box','off','Fontsize',12);
title('In-host Dynamics (Initially with strain 1)','Fontsize',20);
xlabel('Time (years)','Fontsize',16);
ylabel('Frequency','Fontsize',16);
% Legend
legendentries = cell(ns,1);
for i = 1:ns
    legendentries(i,1) = { ['Strain ',num2str(i)] };
end
legend(legendentries,'Location','EastOutside','FontSize',12);

Nfig = subplot(2,2,4);
set(gcf,'DefaultAxesColorOrder',winter(3))
plot(t,N,t,S,t,sum(I),'Linewidth',2.5)
axis([ 0 inf -inf inf]);
set(gca,'box','off','Fontsize',12);
title('','Fontsize',20);
xlabel('Time (years)','Fontsize',16);
ylabel('Population','Fontsize',16);
% Legend
leg = legend('Total', 'Susceptible', 'Infected');
set(leg,'Fontsize',12)



% Total_time = toc(tot);
% disp(['Total time:                                  ',num2str(Total_time,'%.2f'),' s']);
% disp(' ');
