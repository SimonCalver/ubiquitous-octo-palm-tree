%%% FIGURE 3 %%%

% Solves within-host dynamics, and plots the relative frequency of viral
% strains in the active compartment and the latent reservoir.

%close all;
clear all;

%%% PARAMETERS %%%

ns = 2;               % number of strains
s = (1:ns)';           % column vector containing strain indexes
ntsteps = 1e4;         % number of time steps
Tmaxplot = 2;         % maximum time shown in plot

gmin = 1.0;            % replication rate of the least fit strain (per day)
gmax = 1.2;           % replication rate of the fittest strain (per day)
yfact = 365;           % factor to convert from per day to per year
spvlmin = 2;           % minimum log(spvl)
spvlmax = 7;           % maximum log(spvl)

forward_mut = 5e-5;    % Mutation rate of mutating to the next strain (index 1 higher). For neutral evol. use forward_mut = 5e-3
back_mut = 5e-5;       % Mutation rate of mutating into a strain with index 1 lower (NOTE: For non-neutral evolution we set forward_mut = back_mut)
rand_mut = 0;          % Mutation rate of mutating directly to a random other strain (not index 1 higher or lower)

drug = @(t) concentration(t,0.02,0.8,0,1);
g = @(t) yfact * [1.05*(1-drug(t)); 1];
% g = yfact * vir_repr_rate_index_lin(s,gmin,gmax);   % vector with replication rate (per year) for different strains
spvl = spvl_index_lin(s,spvlmin,spvlmax);           % vector with spvl for different strains

%%% Figure setting parameters %%%

plotstrain = 1;        % strain for which the results will be plotted
endfreq = 0.1;         % frequency of initial strain (plotstrain) that is counted as "the timepoint", to measure evolutionary rate

%%% Reservoir parameters %%%

% Probability that newly infected cells become latent, vector for three stages (acute, asymptomatic, AIDS). [Note: in this work we do not change the values of parameters over the three stages]                  
k = 5e-3*[1,1,1];
% Activation rate of latent cells (per day), vector for three stages (acute, asymptomatic, AIDS).
a = 1e-2*[1,1,1];
% Size ratio between latent and active compartment: (size of latent reservoir)/(size of active compartment)
rL = 0.5*[1,1,1];


%%% Duration of infection (as determined by spVL) %%%

D1 = 0.25;                                          % Duration of acute infection (years)
D2 = Hill_func_down(10.^spvl,25.4,3058,0.41);       % Vector of durations of asymptomatic phase (years)
D3 = 0.75;                                          % Duration of AIDS (years)


%%% Time vectors %%%

T1 = D1*ones(ns,1);             % Time of acute infection (vector)
T2 = T1 + D2;                   % Time of acute inf + asymp phase (vector)
T3 = T2 + D3*ones(ns,1);        % Time of complete infection (vector)

Tmax = max(T3)+0.5;             % Maximum time, after this time we stop all calculations
dtau = Tmax/ntsteps;            % time step
tau = 0:dtau:Tmax;              % time vector, with all time steps
ltau = length(tau);             % length of time vector

tau1 = round(ntsteps*T1/Tmax);        % Index of time acute -> asymptomatic
tau2 = round(ntsteps*T2/Tmax);        % Index of time asymptomatic -> AIDS
tau3 = round(ntsteps*T3/Tmax);        % Index of time AIDS -> death


%%% MODEL %%%

tic;

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
for i = plotstrain %1:ns
   [x(:,1:tau1(i),i), y(:,1:tau1(i),i)] = solve_qsODE_PrEP(tau(1:tau1(i)),Init(:,i),ns,g,forward_mut,back_mut,rand_mut,k(1),yfact*a(1),rL(1));       % Acute phase
   disp(['Acute phase solved, strain ',num2str(i,'%d'),'.'])
   tempstart(1:ns) = x(:,tau1(i),i);
   tempstart((ns+1):2*ns) = y(:,tau1(i),i);
   [x(:,(tau1(i)+1):tau2(i),i), y(:,(tau1(i)+1):tau2(i),i)] = solve_qsODE_PrEP(tau((tau1(i)+1):tau2(i)),tempstart,ns,g,forward_mut,back_mut,rand_mut,k(2),yfact*a(2),rL(2));       % Chronic phase
   disp(['Chronic phase solved, strain ',num2str(i,'%d'),'.'])
   tempstart(1:ns) = x(:,tau2(i),i);
   tempstart((ns+1):2*ns) = y(:,tau2(i),i);
   [x(:,tau2(i)+1:end,i), y(:,tau2(i)+1:end,i)] = solve_qsODE_PrEP(tau((tau2(i)+1):end),tempstart,ns,g,forward_mut,back_mut,rand_mut,k(3),yfact*a(3),rL(3));       % AIDS phase
   disp(['AIDS phase solved, strain ',num2str(i,'%d'),'.'])
end

%Calculate at what time the initial strain is reduced to endfreq
endtime = Tmax;
for tt = 1:ltau
    if x(plotstrain,tt,plotstrain) < endfreq
        endtime = tt*dtau;
        break
    end
end

%%% FIGURE %%%


solvetime = toc;
disp(['Time to solve within-host system: ',num2str(solvetime,'%.2f'),' seconds'])

for i=plotstrain %1:ns
    cfig=figure;
    set(cfig,'Units','centimeters','Position',[2 2 35 12])
    clf;
    set(gcf,'DefaultAxesColorOrder',jet(ns))
    
    % Strain distribution in active compartment
    subplot(2,2,1)
    hold on;
    plot(tau(1:end),x(:,1:end,i),'Linewidth',2);
    axis([ 0 Tmaxplot 0 1.05]);
    %semilogy(tau(1:end),x(:,1:end,i),'Linewidth',2);       % to plot on a log y-axis
    %axis([ 0 Tmaxplot 1e-4 1.05]);
%     set(gca,'Xtick',0:2:Tmaxplot,'Fontsize',16);
    title('Active compartment','Fontsize',20);
    xlabel('Time (years)','Fontsize',18);
    ylabel('Strain frequency','Fontsize',18);
    
    
    %Legend
    legendentries = cell(ns,1);
    if gmax == gmin
        %Legend for neutral evolution
        for j = 1:(ns-1)
            legendentries(j,1) = { [num2str(j-1),' mutations'] };
        end
        legendentries(ns,1) = { ['\geq',num2str(ns-1),' mutations'] };
    else
        %Legend for non-neutral evolution
        for j = 1:ns
            legendentries(j,1) = { ['Strain ',num2str(j)] };
        end
    end  
    
    leg = legend(legendentries,'Location','east');
    set(leg,'FontSize',14);
   %}
    
    %Draw a horizontal and vertical line at (endtime, endfreq)
%     line([0 endtime],[endfreq endfreq],'Color','black','Linewidth',1.5);
%     line([endtime endtime],[0 endfreq],'Color','black','Linewidth',1.5);
    
       
    % Strain distribution in reservoir
    subplot(2,2,2)
    hold on;
    plot(tau(1:end),y(:,1:end,i),'Linewidth',2);
    axis([ 0 Tmaxplot 0 1.05]);
    %semilogy(tau(1:end),y(:,1:end,i),'Linewidth',2);       % to plot at a log y axis
    %axis([ 0 Tmaxplot 1e-4 1.05]);
%     set(gca,'Xtick',0:2:Tmaxplot,'Fontsize',16);
    title('Latent reservoir','Fontsize',20);
    xlabel('Time (years)','Fontsize',18);
    ylabel('Strain frequency','Fontsize',18);
    
end


fitness = zeros(ns,length(tau));

for i = 1:length(tau)
    rep = g(tau(i));
    for k = 1:length(s)
    fitness(k,i) = rep(k);
    end
end


ffig=subplot(2,2,4);
set(gcf,'DefaultAxesColorOrder',jet(ns))
plot(tau,fitness,'Linewidth',2.5);
axis([ 0 Tmaxplot -inf inf]);
set(gca,'box','off','Fontsize',12);
title('','Fontsize',20);
xlabel('Time (years)','Fontsize',16);
ylabel('Fitness','Fontsize',16);
% Legend
legendentries = cell(ns,1);
for i = 1:ns
    legendentries(i,1) = { ['Strain ',num2str(i)] };
end
leg = legend(legendentries,'Location','east');
    set(leg,'FontSize',14);


cfig=subplot(2,2,3);
plot(tau,drug(tau),'Linewidth',2.5,'Color',[0.5 0.2 0.8]);
axis([ 0 Tmaxplot -inf inf]);
set(gca,'box','off','Fontsize',12);
title('','Fontsize',20);
xlabel('Time (years)','Fontsize',16);
ylabel('Drug Concentration','Fontsize',16);

