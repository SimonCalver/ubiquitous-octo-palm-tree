%%% Between host dynamics with two host types: on PrEP andnot on PrEP

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
spvlmin = 4;           % minimum log(spvl), assume that these are all the same so we only worry about resisitance to drug, (i.e max = min)?
spvlmax = 4;           % maximum log(spvl)

forward_mut = 5e-5;    % Mutation rate of mutating to the next strain (index 1 higher)
back_mut = 5e-5;       % Mutation rate of mutating into a strain with index 1 lower (NOTE: For non-neutral evolution we set forward_mut = back_mut)
rand_mut = 0;          % Small mutation rate, mutation rate of mutating directly to a random other strain
ART = 2;               % After ART years the host goes onto ART (i.e. is removed from the susceptible population)
PrEP = 0.05;           % proportion of people on PrEP (constant)

drug = @(t) 1;%concentration(t,0.02,0.8,0,1);

% Need two fitness profiles for host on PrEP and not
g = cell(2);

g{1} = @(t) yfact * [1.05; 1]; % Fitness without PrEP
g{2} = @(t) yfact * [1.05*(1-drug(t)); 1]; % Strain fitness within host on PrEP



spvl = spvl_index_lin(s,spvlmin,spvlmax); % Strain set-point viral loads, for now assume they are the samefor each strain

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
D2 = (ART-D1);                                      % Host goes onto ART 

lambda1 = 2.76;                                     % Infectiousness during acute stage
lambda2 = Hill_func_up(10.^spvl,0.317,13938,1.02);  % Vector of infectiousness during asymptomatic phase

%%% WITHIN-HOST DYNAMICS %%%

%%% Time vectors %%%

T1 = D1*ones(ns,1);
T2 = (D1 + D2)*ones(ns,1);

appT1 = round(T1/dt) * dt; % Approximate time of T1 using the coarser grid 
appT2 = round(T2/dt) * dt;
appTmax = max(appT2);
tt = 0:dt:appTmax;
ltt = length(tt);
tau1 = round( appT1 / dt )+1; % index of T1 in the grid
tau2 = round( appT2 / dt )+1; % index of AT2 in the grid
% I have put the "+1"s in because there should be 26 timesteps (when dt =
% 0.01) for D1 but this only gives 25 (I think)

disc_vec = exp(-nu*tt);    % natural mortality over the course of an infection
disc_last = cell(2); % The last point of the vector of discount factor due to natural mortality
disc_last{1} = exp(-nu*appT2);% At the moment we assume that all infected are removed (cured) after a fixed time so these values are the same for each host type.
disc_last{2} = exp(-nu*appT2);% This need not be the case

host_Alphas = cell(2);
for i = 1:2
    host_Alphas{i} = zeros(ns,ltt);
end

for j = 1:2
    for i = 1:ns
        host_Alphas{j}(i,1:tau1(i)) = lambda1;      % I chose that the step function is right continuous (the value in an interval is that at the right end). The acute phase has no area, if lT1=1
        host_Alphas{j}(i,(tau1(i)+1):tau2(i)) = lambda2(i);
        host_Alphas{j}(i,:) = Host_Infectivity(j,i,host_Alphas{j}(i,:));
        host_Alphas{j}(i,:) = host_Alphas{j}(i,:) .* disc_vec;      % Incorporate natural mortality in total infectivity, because these terms always occur together
    end
end


% Initial conditions - start with 1 infecting virion (in each compartment)
Init_temp = eye(ns);
Init = zeros(2*ns,ns);
Init(1:ns, :) = Init_temp;
Init((ns+1):2*ns, :) = Init_temp;

% Solve within-host system, using ODE solver

host_x = cell(2);
host_y = cell(2);

for i = 1:2
    host_x{i} = zeros(ns,ltt,ns); % One 2D matrix for each starting point
    host_y{i} = zeros(ns,ltt,ns); % One 2D matrix for each starting point
end

tempstart = zeros(2*ns,1);

for j = 1:2
    for i = 1:ns
       [host_x{j}(:,1:tau1(i),i), host_y{j}(:,1:tau1(i),i)] = solve_qsODE_PrEP(tt(1:tau1(i)),Init(:,i),ns,g{j},forward_mut,back_mut,rand_mut,k(1),yfact*a(1),rL(1));       % Acute phase
       tempstart(1:ns) = host_x{j}(:,tau1(i),i);
       tempstart((ns+1):2*ns) = host_y{j}(:,tau1(i),i);
       [host_x{j}(:,(tau1(i)+1):tau2(i),i), host_y{j}(:,(tau1(i)+1):tau2(i),i)] = solve_qsODE_PrEP(tt((tau1(i)+1):tau2(i)),tempstart,ns,g{j},forward_mut,back_mut,rand_mut,k(2),yfact*a(2),rL(2));       % Chronic phase
     end
end


% Determine strain-specific infectivity profiles 
% betas = zeros(ns,ltt,ns);   % strain-specific infectivity profiles (with natural mortality incorporated)
% Kc = zeros(ns);             % Next generation matrix
% for i = 1:ns
%     for j = 1:ns
%         for itt = 1:tau3(j)
%             betas(i,itt,j) = alphas(j,itt) * x(i,itt,j) ; 
%         end
%         Kc(i,j) = dt * trapz( betas(i,1:tau3(j),j) );
%     end
% end

%%% BETWEEN-HOST DYNAMICS %%%

%%% Initialisation %%%

ICint = appTmax; % Duration of interval of initial conditions, before tmin, during which incidence is constant at I0
iT0 = ICint/dt;

t = (tmin-ICint):dt:tmax;
lt = length(t);
% I = zeros(ns,lt);       % Matrix of incidence over absolute time, one row per type
% Ibar = zeros(ns,lt);    % Matrix of relative incidence over time
% Itot = zeros(1,lt);     % Total incidence
% J = zeros(ns,lt);       % Prevalence is 0 at the beginning, because there are only newly infected initial infective
% Jbar = zeros(ns,lt);    % Relative prevalence
% Jtot = zeros(1,lt);     % Total prevalence
% S = zeros(1,lt);        % Total susceptibles are N at the beginning, for the same reason
% Y = zeros(ns,lt);       % Total amount of virus in the population, one row per strain
% N = zeros(1,lt);        % Total population size, which is variable and starts at N0

setup_1 = cell(2);
setup_ns = cell(2);
for i = 1:2
    setup_1{i} = zeros(1,lt);
    setup_ns{i} = zeros(ns,lt);
end

I = setup_ns;       % Matrix of incidence over absolute time, one row per type
Ibar = setup_ns;    % Matrix of relative incidence over time
Itot = setup_1;     % Total incidence
J = setup_ns;       % Prevalence is 0 at the beginning, because there are only newly infected initial infective
Jbar = setup_ns;    % Relative prevalence
Jtot = setup_1;     % Total prevalence
S = setup_1;        % Total susceptibles are N at the beginning, for the same reason
Y = setup_ns;       % Total amount of virus in the population, one row per strain
N = zeros(1,lt);    % Total population size, which is variable and starts at N0

%%%initial conditions %%%


% Basic initial conditions
I0 = cell(2);
for i = 1:2
    I0{i} = zeros(ns,1);   % Column vector of initial values for each strain
end

I0{1}(startstrain) = 1;          % Initial condition for the incidence in the population, with only strain 1 present, we assume it starts in someone not on PrEP
S0 = N0;

for iic = 1:iT0
    for st = 1:2
        I{st}(:,iic) = I0{st};
        Ibar{st}(:,iic) = I0{st} / (sum(I0{1}) + sum(I0{2})); % maybe don't need totals in different hosts
        Itot{st}(iic) = sum(I0{st});
        J{st}(:,iic) = zeros(ns,1);         % Could do better, probably
        Jbar{st}(:,iic) = zeros(ns,1); 
        Jtot{st}(iic) = 0; 
    end
    N(iic) = N0;
    S{1}(iic) = (1-PrEP)*S0; % susceptiples not on PrEP 
    S{2}(iic) = PrEP*S0;
end

N(iT0) = N0;
S{1}(iT0) = (1-PrEP) * (N(iT0) - (Itot{1}(iT0)+Itot{2}(iT0)) * dt * trapz( disc_vec ));
S{2}(iT0) = PrEP * (N(iT0) - (Itot{1}(iT0)+Itot{2}(iT0)) * dt * trapz( disc_vec ));


% Solve between-host system
FOI = cell(2);         % Force of infection of different strains over time for different hosts
FOI{1} = zeros(ns,lt);
FOI{2} = zeros(ns,lt);
FOI_temp2 = cell(2);
FOI_temp2{1} = zeros(ns,1);
FOI_temp2{2} = zeros(ns,1);
for it = (iT0+1):lt
    discI = cell(2);
    discI{1} = zeros(ns,1);
    discI{2} = zeros(ns,1);
    for i = 1:ns
        FOI_temp = cell(2);
        FOI_temp{1} = zeros(ns,ltt);       % FOI from type j at any point in the past (acting on i, but not depending explicitly, i.e. overwritten every loop)
        FOI_temp{2} = zeros(ns,ltt); 
        for j = 1:ns
            FOI_temp{1}(j,1) = 0;          % The first point of this vector is always 0, because it refers to the time point t(it), the one we are calculating the solution at
            FOI_temp{2}(j,1) = 0;
            for itt = 2:min(it,tau2(j))
                FOI_temp{1}(j,itt) = host_Alphas{1}(j,itt) * host_x{1}(i,itt,j) * I{1}(j,it-itt+1) +...
                                     host_Alphas{1}(j,itt) * host_x{2}(i,itt,j) * I{2}(j,it-itt+1); % at the moment we are integrating over the same time, so these can go together
                FOI_temp{2}(j,itt) = host_Alphas{2}(j,itt) * host_x{1}(i,itt,j) * I{1}(j,it-itt+1) +...
                                     host_Alphas{2}(j,itt) * host_x{2}(i,itt,j) * I{2}(j,it-itt+1); % at the moment we are integrating over the same time, so these can go together
            end
            FOI_temp2{1}(j) = dt * trapz( FOI_temp{1}(j,1:min(it,tau2(j))) );     % Contibution of j to the FOI acting on i
            FOI_temp2{2}(j) = dt * trapz( FOI_temp{2}(j,1:min(it,tau2(j))) );
        end
        FOI{1}(i,it) = sum( FOI_temp2{1} );
        FOI{2}(i,it) = sum( FOI_temp2{2} );
        I{1}(i,it) = S{1}(it-1) * FOI{1}(i,it) / N(it-1);
        I{2}(i,it) = S{2}(it-1) * FOI{2}(i,it) / N(it-1);
        J{1}(i,it) = dt * trapz( I{1}(i,max(1,it:-1:it-tau2(i)+1)) .* disc_vec(1:tau2(i)) );
        J{2}(i,it) = dt * trapz( I{2}(i,max(1,it:-1:it-tau2(i)+1)) .* disc_vec(1:tau2(i)) );
        discI{1}(i) = I{1}(i,it-tau2(i)+1) * disc_last{1}(i);
        discI{1}(i) = I{1}(i,it-tau2(i)+1) * disc_last{1}(i);% beause of the "+1"s I added earlier I need to take 1 from tau3 to make it the same. 
        % This may well be wrong. These lines used to look like this:
        % J(i,it) = dt * trapz( I(i,max(1,it:-1:it-tau3(i))) .* disc_vec(1:tau3(i)+1) ); % I suspect the max is not necessary
        % discI(i) = I(i,it-tau3(i)) * disc_last(i);
    end
    Itot{1}(it) = sum(I{1}(:,it));
    Itot{2}(it) = sum(I{2}(:,it));
    Jtot{1}(it) = sum(J{1}(:,it));
    Jtot{2}(it) = sum(J{2}(:,it));
    S = N(it-1) - Jtot(it);
    S{1}(it) = (1-PrEP) * S;
    S{2}(it) = (1-PrEP) * S;
    N(it) = N(it-1) + dt * ( B - nu * N(it-1) - ( sum( discI{1} ) + sum( discI{2} ) ) );
    Ibar{1}(:,it) = I{1}(:,it) / ( Itot{1}(it) +  Itot{2}(it) );
    Ibar{2}(:,it) = I{2}(:,it) / ( Itot{1}(it) +  Itot{2}(it) );
    Jbar{1}(:,it) = J{1}(:,it) / ( Jtot{1}(it) + Jtot{2}(it) );
    Jbar{2}(:,it) = J(:,it) / ( Jtot{1}(it) + Jtot{2}(it) );
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
