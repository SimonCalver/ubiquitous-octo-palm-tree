%%% Between host dynamics with two host types: on PrEP andnot on PrEP

% Calculate the strain distribution in a host population

% close all;
clear all;

%%% PARAMETERS %%%

ns = 2;             % number of strains
nh = 2;             % number of host types
s = (1:ns)';        % Column vector containing the strain index
startstrain = 1;    % Strain from which the infection is started
starthost = 2;      % Host from which the infection is started

yfact = 365;           % factor to convert from per day to per year
spvlmin = 5;           % minimum log(spvl), assume that these are all the same so we only worry about resisitance to drug, (i.e max = min)?
spvlmax = 4.5;         % maximum log(spvl)

forward_mut = 5e-5;    % Mutation rate of mutating to the next strain (index 1 higher)
back_mut = 5e-5;       % Mutation rate of mutating into a strain with index 1 lower (NOTE: For non-neutral evolution we set forward_mut = back_mut)
rand_mut = 0;          % Small mutation rate, mutation rate of mutating directly to a random other strain
ART = 2;               % After ART years the host goes onto ART (i.e. is removed from the susceptible population)

Hprops = [0.0 1.0];% 0.0];   % proportions of hosts of each type (must add to 1)
Drug_Adherance = {@(t) 0 @(t) 1.0};% @(t) concentration(t,7/365,0.5,48/8760,0,1.0)};  %concentration(t,0.02,0.8,0,1);
relInfectivity = [1.0, 0.01, 0.3;
                  1.0, 1.0,  1.0];
                  
g = cell(nh,1);

for i = 1:nh
    g{i} = @(t) yfact * [1.05*(1-Drug_Adherance{i}(t)); 1.0];
end

spvl = spvl_index_lin(s,spvlmin,spvlmax); % Strain set-point viral loads, for now assume they are the samefor each strain

N0 = 10000; % Total initial population size
B = 200; % Constant total birth rate per year
nu = 0.02; % Constant per capita yearly death rate

dt = 0.01; % Time-step for the computation of the full dynamics
tmin = 0;
tmax = 300;


%%% Reservoir parameters %%%

% Probability that newly infected cells become latent, vector for three stages (acute, asymptomatic, AIDS). [Note: in this work we do not change the values of parameters over the three stages]                  
k = 0.0051*[1,1,1];
% Activation rate of latent cells (per day), vector for three stages (acute, asymptomatic, AIDS).
a = 1e-2*[1,1,1];
% Size ratio between latent and active compartment: (size of latent reservoir)/(size of active compartment)
rL = 2*[1,1,1];


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
disc_last = cell(nh,1); % The last point of the vector of discount factor due to natural mortality

host_Alphas = cell(nh,1);
for i = 1:nh
    host_Alphas{i} = zeros(ns,ltt);
    disc_last{i} = exp(-nu*appT2); % At the moment we assume that all infected are removed (cured) after a fixed time so these values are the same for each host type.
                                   % This need not be the case
end

for j = 1:nh
    for i = 1:ns
        host_Alphas{j}(i,1:tau1(i)) = lambda1;      % I chose that the step function is right continuous (the value in an interval is that at the right end). The acute phase has no area, if lT1=1
        host_Alphas{j}(i,(tau1(i)+1):tau2(i)) = lambda2(i);
        host_Alphas{j}(i,:) = host_Alphas{j}(i,:) * relInfectivity(i,j); %Host_Infectivity(j,i,host_Alphas{j}(i,:));
        host_Alphas{j}(i,:) = host_Alphas{j}(i,:) .* disc_vec;      % Incorporate natural mortality in total infectivity, because these terms always occur together
    end
end



% Initial conditions - start with 1 infecting virion (in each compartment)
Init_temp = eye(ns);
Init = zeros(2*ns,ns);
Init(1:ns, :) = Init_temp;
Init((ns+1):2*ns, :) = Init_temp;

% Solve within-host system, using ODE solver

host_x = cell(nh,1);
host_y = cell(nh,1);

for i = 1:nh
    host_x{i} = zeros(ns,ltt,ns); % One 2D matrix for each starting point
    host_y{i} = zeros(ns,ltt,ns); % One 2D matrix for each starting point
end

tempstart = zeros(2*ns,1);

for j = 1:nh
    for i = 1:ns
       [host_x{j}(:,1:tau1(i),i), host_y{j}(:,1:tau1(i),i)] = solve_qsODE_PrEP(tt(1:tau1(i)),Init(:,i),ns,g{j},forward_mut,back_mut,rand_mut,k(1),yfact*a(1),rL(1));       % Acute phase
       tempstart(1:ns) = host_x{j}(:,tau1(i),i);
       tempstart((ns+1):2*ns) = host_y{j}(:,tau1(i),i);
       [host_x{j}(:,(tau1(i)+1):tau2(i),i), host_y{j}(:,(tau1(i)+1):tau2(i),i)] = solve_qsODE_PrEP(tt((tau1(i)+1):tau2(i)),tempstart,ns,g{j},forward_mut,back_mut,rand_mut,k(2),yfact*a(2),rL(2));       % Chronic phase
     end
end


%%%%%%%%%% Viral Load
%%%%%%%%%%%%%%%%%%%%%

VL = cell(nh,1);
for i = 1:nh
    VL{i} = zeros(ns,length(tt),ns);
end

for j = 1:nh
    for i = 1:ns
        for k = 1:ns
        VL{j}(k,1:tau1(i),i) = lambda1 * host_x{j}(k,1:tau1(i),i);
        VL{j}(k,(tau1(i)+1):tau2(i),i) = lambda2(i) * host_x{j}(k,(tau1(i)+1):tau2(i),i);
        end
    end
end    

% Strain 1 VL is reduced by the drug, so it is ~0 when a lot of drug, when
% initiated by either strain. It won't go completely to 0 (andif it did
% there would be NaNs in the solution). When infected initiLLY WITH STRAIN
% 2, the VL follows the dynnamics in the absence of drug
for i = 1:nh
    VL{i}(1,1:end,1) =  max(VL{i}(1,1:end,1) .* (1 - Drug_Adherance{i}(tt)), 0.0001);
end
% Lets say after 30 days the VL of the resistant strain reaches the spvl,
% since there is definitely some present 
tStart = find(Drug_Adherance{i}(tt) >= 0, 1);
tStop = find(tt >= tt(tStart) + 0.0822, 1);

for i = 1:nh
    VL{i}(2,tStart:tStop,1) =  min(VL{i}(2,tStart:tStop,1) .* (1 - Drug_Adherance{i}(tt(tStart:tStop))), 0.0001);
end

% Incorporate the natural mortality
for j = 1:nh
    for i = 1:2
         for k = 1:2
            VL{j}(k,1:end,i) = VL{j}(k,1:end,i) .* exp(-nu*tt);     
         end
    end
end

% VL{2}(1,1:end,1) = 0.00001 * VL{2}(1,1:end,1);
% VL{2}(2,1:end,1) = 0.00001 * VL{2}(2,1:end,1);
% VL{2}(1,1:end,2) = 0.00001 * VL{2}(1,1:end,2);
% VL{2}(2,1:end,2) = 0.00001 * VL{2}(2,1:end,2);

% In VL{k}(i,t,j) strain i inititaes infection in new host, so this is
% infctivity of strain i in host k when donor was initially infected with 
% strain j

%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%



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

setup_1 = cell(nh,1);
setup_ns = cell(nh,1);
for i = 1:nh
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
I0 = cell(nh,1);
for i = 1:nh
    I0{i} = zeros(ns,1);   % Column vector of initial values for each strain
end

I0{starthost}(startstrain) = 1;          % Initial condition for the incidence in the population, with only strain 1 present, we assume it starts in someone not on PrEP
S0 = N0;

for iic = 1:iT0
    for st = 1:nh
        I{st}(:,iic) = I0{st};
        Isum = 0;%zeros(iT0,1);
        for i = 1:nh 
            Isum = Isum + sum(I0{i});
        end
        Ibar{st}(:,iic) = I0{st} / Isum; % maybe don't need totals in different hosts
        Itot{st}(iic) = sum(I0{st});
        J{st}(:,iic) = zeros(ns,1);         % Could do better, probably
        Jbar{st}(:,iic) = zeros(ns,1); 
        Jtot{st}(iic) = 0; 
        S{st}(iic) = Hprops(st)*S0; 
    end
    N(iic) = N0;
end

N(iT0) = N0;
for i = 1:nh
    S{i}(iT0) = Hprops(i) * (N(iT0) - Isum * dt * trapz( disc_vec ));
end
%     
% S{1}(iT0) = Hprops(1) * (N(iT0) - (Itot{1}(iT0)+Itot{2}(iT0)) * dt * trapz( disc_vec ));
% S{2}(iT0) = Hprops(2) * (N(iT0) - (Itot{1}(iT0)+Itot{2}(iT0)) * dt * trapz( disc_vec ));


% Solve between-host system
FOI = cell(nh,1);         % Force of infection of different strains over time for different hosts
FOI_temp2 = cell(nh,1);
for i = 1:nh
    FOI{i} = zeros(ns,lt);
    FOI_temp2{i} = zeros(ns,1);
end


for it = (iT0+1):lt
    discI = cell(nh,1);
    for discI_ind = 1:nh
        discI{discI_ind} = zeros(ns,1);
    end
    for i = 1:ns % The strain that initiates infection in new host 
        FOI_temp = cell(nh,1);  % FOI from type j at any point in the past (acting on i, but not depending explicitly, i.e. overwritten every loop)
        for FOI_temp_ind = 1:nh
            FOI_temp{FOI_temp_ind} = zeros(ns,ltt);
        end       
        for j = 1:ns % The strain that initiated infection in the donor (determines VL)
            for FOI_temp_ind = 1:nh
                FOI_temp{FOI_temp_ind}(j,1) = 0;    % The first point of this vector is always 0, because it refers to the time point t(it), the one we are calculating the solution at
            end
            for itt = 2:min(it,tau2(j))
                for FOI_temp_ind1 = 1:nh
                    for FOI_temp_ind2 = 1:nh
                        % The force of infection of virus initiated by strain j for each host
                        FOI_temp{FOI_temp_ind1}(j,itt) = FOI_temp{FOI_temp_ind1}(j,itt) + VL{FOI_temp_ind1}(i,itt,j) * I{FOI_temp_ind2}(j,it-itt+1);% host_x{FOI_temp_ind2}(i,itt,j) * I{FOI_temp_ind2}(j,it-itt+1);
                    end
%                     FOI_temp{FOI_temp_ind1}(j,itt) = FOI_temp{FOI_temp_ind1}(j,itt) * host_Alphas{FOI_temp_ind1}(j,itt);      
                end
            end
                
            for FOI_temp2_ind = 1:nh
                FOI_temp2{FOI_temp2_ind}(j) = dt * trapz( FOI_temp{FOI_temp2_ind}(j,1:min(it,tau2(j))) );     % Contibution of j to the FOI acting on i
            end
        end
        for host = 1:nh
            FOI{host}(i,it) = sum( FOI_temp2{host} );%   FOI_temp2{host}(i);%      should this be a sum? It was already summed before the integral
            I{host}(i,it) = S{host}(it-1) * FOI{host}(i,it) / N(it-1);
            J{host}(i,it) = dt * trapz( I{host}(i,max(1,it:-1:it-tau2(i)+1)) .* disc_vec(1:tau2(i)) );
            discI{host}(i) = I{host}(i,it-tau2(i)+1) * disc_last{host}(i);
        end      
%         FOI{1}(i,it) = sum( FOI_temp2{1} );
%         FOI{2}(i,it) = sum( FOI_temp2{2} );
%         I{1}(i,it) = S{1}(it-1) * FOI{1}(i,it) / N(it-1);
%         I{2}(i,it) = S{2}(it-1) * FOI{2}(i,it) / N(it-1);
%         J{1}(i,it) = dt * trapz( I{1}(i,max(1,it:-1:it-tau2(i)+1)) .* disc_vec(1:tau2(i)) );
%         J{2}(i,it) = dt * trapz( I{2}(i,max(1,it:-1:it-tau2(i)+1)) .* disc_vec(1:tau2(i)) );
%         discI{1}(i) = I{1}(i,it-tau2(i)+1) * disc_last{1}(i);
%         discI{2}(i) = I{2}(i,it-tau2(i)+1) * disc_last{2}(i);% beause of the "+1"s I added earlier I need to take 1 from tau3 to make it the same. 
        % This may well be wrong. These lines used to look like this:
        % J(i,it) = dt * trapz( I(i,max(1,it:-1:it-tau3(i))) .* disc_vec(1:tau3(i)+1) ); % I suspect the max is not necessary
        % discI(i) = I(i,it-tau3(i)) * disc_last(i);
    end
    Isum = 0;
    Jsum = 0;
    discIsum = 0;
    for host = 1:nh
        Itot{host}(it) = sum(I{host}(:,it));
        Jtot{host}(it) = sum(J{host}(:,it));
    end
    for ind = 1:nh 
        Isum = Isum + sum(Itot{ind}(it));
        Jsum = Jsum + sum(Jtot{ind}(it));
        discIsum = discIsum + sum(discI{ind});
    end
    
    
    if(it == 10000)
        Hprops = [1.0 0.0 0.0]; % introduce PrEP
    end
    
    Sus = N(it-1) - Jsum;
    for host = 1:nh
%%%%%%%%% For fixed proportion use this
        S{host}(it) = Hprops(host) * Sus;
%%%%%%%%%
%%%%%%%%% for a fixed initial amount use this THIS IS WRONG
%         S{host}(it) = S{host}(it-1) + dt * ( Hprops(host) * B - nu * S{host}(it-1) - sum(discI{host}) );
%%%%%%%%% 
        N(it) = N(it-1) + dt * ( B - nu * N(it-1) - discIsum );
        Ibar{host}(:,it) = I{host}(:,it) / Isum;
        Jbar{host}(:,it) = J{host}(:,it) / Jsum;
    end
    
%%%%% This is stupid
    N(N < 1) = 0e-20;%  N = floor(N);
    S{1}(S{1} < 1) = 0e-20;
    S{2}(S{2} < 1) = 0e-20;%  S{2} = floor(S{2});
    J{1}(I{1} < 1) = 0e-20;%  I{1} = floor(I{1});
    J{2}(I{2} < 1) = 0e-20;%  I{2} = floor(I{2});
%%%%%
  
end

Prevalence = Jbar{1};
Susceptibles = S{1};
Infected = sum(J{1});
for i = 2:nh
    Prevalence = Prevalence + Jbar{i};
    Susceptibles = Susceptibles + S{i};
    Infected = Infected + sum(J{i});
end

%% FIGURES 


% scrsz = get(groot,'ScreenSize');
% main = figure('Position',[1 0 scrsz(3)/2 scrsz(4)/2]);

figure

% cfig = subplot(1,2,1);
set(gca,'DefaultAxesColorOrder',winter(ns))

plot(t(find(t==0):end),Prevalence(:,find(t==0):end),'Linewidth',2.5); % The find makes the figure look nicer, but maybe it would be better to change Jbar so it isn't zero for t<0?
axis([ 0 tmax 0 1]);
set(gca,'box','off','Fontsize',18);
% title('Stratified prevalence','Fontsize',22,'interpreter','latex');
xlabel('Time (years)','Fontsize',20,'interpreter','latex');
ylabel('Strain frequency','Fontsize',20,'interpreter','latex');
legend({'Wild-type' 'Resistant'},'Location','East','FontSize',20,'interpreter','latex');



% Nfig = subplot(1,2,2);
figure
set(gca,'DefaultAxesColorOrder',summer(nh+1))
hold on;
% legendentries = cell(nh+1,1);
% for i = 1:nh
%     plot(t,S{i},'Linewidth',2.5);
%     legendentries(i,1) = { ['Susceptible ',num2str(i)] };
% end
% legendentries(nh+1,1) = { 'Infected ' };
plot(t,Susceptibles,'Linewidth',2.5); hold on;
legendentries = {'Susceptible', 'Infected'};
plot(t,Infected,'Linewidth',2.5); 
axis([ 0 tmax 0 inf]);
set(gca,'box','off','Fontsize',18);
% title('Population','Fontsize',22,'interpreter','latex');
xlabel('Time (years)','Fontsize',20,'interpreter','latex');
ylabel('Population','Fontsize',20,'interpreter','latex');
legend(legendentries,'Location','East','FontSize',20,'interpreter','latex');


%%% FIGURES %%%
% 
% 
% scrsz = get(groot,'ScreenSize');
% main = figure('Position',[1 0 scrsz(3) scrsz(4)]);
% 
% 
% 
% cfig = subplot(2,2,1);
% set(gcf,'DefaultAxesColorOrder',jet(ns))
% plot(t,Jbar,'Linewidth',2.5);
% axis([ 0 tmax 0 1.005]);
% set(gca,'XTick',0:25:tmax,'YTick',0:0.2:1,'box','off','Fontsize',12);
% title('Stratified prevalence','Fontsize',20);
% xlabel('Time (years)','Fontsize',16);
% ylabel('Strain frequency','Fontsize',16);
% % Legend
% legendentries = cell(ns,1);
% for i = 1:ns
%     legendentries(i,1) = { ['Strain ',num2str(i)] };
% end
% legend(legendentries,'Location','EastOutside','FontSize',12);
% str(1) = {['R0 = ',num2str(R0_end)]};
% str(2) = {['log(spVL) = ',num2str(avSPVL_end)]};
% R0text = text(80,0.6,str,'FontSize',14);
% set(R0text,'EdgeColor',[0 0 0])
% 
% 
% fitness = zeros(ns,length(tt));
% 
% for i = 1:length(tt)
%     rep = g(tt(i));
%     for k = 1:length(s)
%     fitness(k,i) = rep(k);
%     end
% end
% 
% 
% ffig=subplot(2,2,2);
% set(gcf,'DefaultAxesColorOrder',jet(ns))
% plot(tt,fitness,'Linewidth',2.5);
% axis([ 0 inf -inf inf]);
% set(gca,'box','off','Fontsize',12);
% title('','Fontsize',20);
% xlabel('Time (years)','Fontsize',16);
% ylabel('Fitness','Fontsize',16);
% % Legend
% legendentries = cell(ns,1);
% for i = 1:ns
%     legendentries(i,1) = { ['Strain ',num2str(i)] };
% end
% legend(legendentries,'Location','EastOutside','FontSize',12);
% 
% 
% 
% ifig = subplot(2,2,3);
% set(gcf,'DefaultAxesColorOrder',jet(ns))
% plot(tt,x(:,:,1),'Linewidth',2.5);
% axis([ 0 inf -inf inf]);
% set(gca,'box','off','Fontsize',12);
% title('In-host Dynamics (Initially with strain 1)','Fontsize',20);
% xlabel('Time (years)','Fontsize',16);
% ylabel('Frequency','Fontsize',16);
% % Legend
% legendentries = cell(ns,1);
% for i = 1:ns
%     legendentries(i,1) = { ['Strain ',num2str(i)] };
% end
% legend(legendentries,'Location','EastOutside','FontSize',12);
% 
% Nfig = subplot(2,2,4);
% set(gcf,'DefaultAxesColorOrder',winter(3))
% plot(t,N,t,S,t,sum(I),'Linewidth',2.5)
% axis([ 0 inf -inf inf]);
% set(gca,'box','off','Fontsize',12);
% title('','Fontsize',20);
% xlabel('Time (years)','Fontsize',16);
% ylabel('Population','Fontsize',16);
% % Legend
% leg = legend('Total', 'Susceptible', 'Infected');
% set(leg,'Fontsize',12)
% 
% load handel;
% sound(y,Fs); % Check this out http://ssli.ee.washington.edu/courses/ee299/labs/Lab3.pdf

