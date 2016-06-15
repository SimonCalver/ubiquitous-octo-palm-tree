%%% Between host dynamics with two host types: on PrEP andnot on PrEP

% Calculate the strain distribution in a host population

% close all;
clear all;

%%% PARAMETERS %%%

ns = 2;             % number of strains
nh = 4;             % number of host types
s = (1:ns)';        % Column vector containing the strain index
startstrain = 1;    % Strain from which the infection is started
starthost = 1;      % Host from which the infection is started

yfact = 365;           % factor to convert from per day to per year
spvlmin = 7;           % minimum log(spvl), assume that these are all the same so we only worry about resisitance to drug, (i.e max = min)?
spvlmax = 7;           % maximum log(spvl)

forward_mut = (5e-5)^1;    % Mutation rate of mutating to the next strain (index 1 higher)
back_mut = (5e-5)^1;       % Mutation rate of mutating into a strain with index 1 lower (NOTE: For non-neutral evolution we set forward_mut = back_mut)
rand_mut = 0;              % Small mutation rate, mutation rate of mutating directly to a random other strain

inf_test = 0.5;        % Time after infection that virus is detected for hosts on PrEP
ART_PrEP = 1.0;        % Time after infection is discovered that host goes onto ART (i.e. is removed from the susceptible population)
ART = 2.0;             % The time hosts not on PrEP get treatment 
Delay = 20/365;        % Time after infection that host goes onto PrEP

% Hprobs = [0.0 1.0 0.0 0.0];
Hprobs = [0.8 0.1560 0.004 0.04];   % proportions of hosts of each type (must add to 1)
Hprobs2 = Hprobs;

transmission = 0.0; % Proportion of virus that is transmitted from the reservoir

duh = 1.2;

Drug_Adherance = {@(t) 0*t  @(t) concentration(t,0,inf_test,48/8760,0,1.0)  @(t) concentration(t,Delay,inf_test,48/8760,0,1.0) @(t) concentration(t,0,inf_test,48/8760,0,0.8)};  %concentration(t,0.02,0.8,0,1);
  
% Here's a thought; what if the amount of drug depends on wether or not the person is infected? so one determines their infectvity and the other how easily they can be infected                    
Drug_Adherance_Uninfected = {@(t) 0*t  @(t) 1.0+0*t  @(t) 1.0*(t >= Delay) @(t) 0.8+0*t};  %concentration(t,0.02,0.8,0,1);


g = cell(nh,1);

for i = 1:nh
    g{i} = @(t) yfact * [1.0*(1-Drug_Adherance{i}(t)); 0.75];
end

spvl = spvl_index_lin(s,spvlmin,spvlmax); % Strain set-point viral loads, for now assume they are the samefor each strain

N0 = 10000; % Total initial population size
B = 200; % Constant total birth rate per year
nu = 0.02; % Constant per capita yearly death rate

dt = 0.01; % Time-step for the computation of the full dynamics
tmin = 0;
tmax = 100;


%%% Reservoir parameters %%%

% Probability that newly infected cells become latent, vector for three stages (acute, asymptomatic, AIDS). [Note: in this work we do not change the values of parameters over the three stages]                  
k = 0.0052*[1,1,1];%  0.016*[1,1,1];%   0.0051
% Activation rate of latent cells (per day), vector for three stages (acute, asymptomatic, AIDS).
a = 1e-2*[1,1,1];%    0.003*[1,1,1];% 1e-2  
% Size ratio between latent and active compartment: (size of latent reservoir)/(size of active compartment)
rL = 2*[1,1,1];%  0.95*[1,1,1];%   2


%% Predefined infectivity profile

% Multiply these up to make R0 larger
lambda1 = duh * 2.76;                                     % Infectiousness during acute stage
lambda2 = duh * Hill_func_up(10.^spvl,0.317,13938,1.02);  % Vector of infectiousness during asymptomatic phase

%%% WITHIN-HOST DYNAMICS %%%

%%% Time vectors %%%

T1 = {0.25*ones(ns,1) 0.25*ones(ns,1) 0.25*ones(ns,1) 0.25*ones(ns,1)};
T2 = {ART*ones(ns,1) (inf_test+ART_PrEP)*ones(ns,1) (inf_test+ART_PrEP+Delay)*ones(ns,1) (inf_test+ART_PrEP)*ones(ns,1)};

    appT1 = cell(nh,1);  
    appT2 = cell(nh,1); 
    
for i1 = 1:nh
    appT1{i1} = round(T1{i1}/dt) * dt; % Approximate time of T1 using the coarser grid 
    appT2{i1} = round(T2{i1}/dt) * dt;
end

appTmax = 0;
for i1 = 1:nh
    if appTmax < max(appT2{i1})
        appTmax = max(appT2{i1});
    end
end

% % appTmax = max(appT2);
tt = 0:dt:appTmax;
ltt = length(tt);

    tau1 = cell(nh,1);  
    tau2 = cell(nh,1); 
    
for i1 = 1:nh
    tau1{i1} = round( appT1{i1} / dt ) + 1; % index of T1 in the grid
    tau2{i1} = round( appT2{i1} / dt ) + 1; % index of AT2 in the grid
% I have put the "+1"s in because there should be 26 timesteps (when dt =
% 0.01) for D1 but this only gives 25 (I think)
end

disc_vec = exp(-nu*tt);    % natural mortality over the course of an infection
disc_last = cell(nh,1); % The last point of the vector of discount factor due to natural mortality




% Initial conditions - start with 1 infecting virion (in each compartment)
Init_temp = eye(ns);
Init = zeros(2*ns,ns);
Init(1:ns, :) = Init_temp;
Init((ns+1):2*ns, :) = Init_temp;

% Ignore this!!
% InitTest = 0.5*ones(4,2);
% Init = InitTest;

% Solve within-host system, using ODE solver

host_x = cell(nh,1);
host_y = cell(nh,1);

for i = 1:nh
    host_x{i} = zeros(ns,ltt,ns); % One 2D matrix for each starting point
    host_y{i} = zeros(ns,ltt,ns); % One 2D matrix for each starting point
    disc_last{i} = exp(-nu*appT2{i});
end

for j = 1:nh
    for i = 1:ns
       [host_x{j}(:,1:tau2{j}(i),i), host_y{j}(:,1:tau2{j}(i),i)] = solve_qsODE_PrEP(tt(1:tau2{j}(i)),Init(:,i),ns,g{j},forward_mut,back_mut,rand_mut,k(1),yfact*a(1),rL(1));
    end
end


%% Infectivity
%%%%%%%%%%%%%%%%%%%%

% From the frequencies find frequency in each host, this will depend on
% drug concentration and stuff

c0 = 0.9; % Competition parameter without the drug, when this is 1 WT
% cannot cause infection when there is no drug
c1 = 0.5; % Competition parameter with the drug, when this is 1 WT 
% (susceptible) cannot cause infection in host on PrEP

x_adjusted = cell(nh,nh);
y_adjusted = cell(nh,nh);
for i = 1:nh
    for j = 1:nh
        x_adjusted{i,j} = zeros(ns,length(tt),ns);
        y_adjusted{i,j} = zeros(ns,length(tt),ns);
    end
end
% Simple way of allowing for competition, best if these don't sum to one.
% So for recipient on the druginfectivity of each is lowered compared to
% PrEP naive host. Transmissibility of both is affectedby the drug either
% positively or negatively
for i1 = 1:nh
    for i2 = 1:nh
        for j1 = 1:ns
             for j2 = 1:ns
                x_adjusted{i1,i2}(j2,1:tau2{i2}(j1),j1) = ( (1-transmission) * host_x{i2}(j2,1:tau2{i2}(j1),j1) + transmission * host_y{i2}(j2,1:tau2{i2}(j1),j1) );
             end
        end
    end
end

          
% cR = 0.1; % Vary from about 0.1 to infinity
% cWT = 0.5;
% 
% for i1 = 1:nh
%     for i2 = 1:nh
%         for j1 = 1:ns             
%             c = Drug_Adherance_Uninfected{i1}(tt);%c0*(1 - Drug_Adherance_Uninfected{i1}(tt)) + c1*Drug_Adherance_Uninfected{i1}(tt);
% 
% %             x_adjusted{i1,i2}(1,1:end,j1) = 1/50*(x_adjusted{i1,i2}(1,1:end,j1)+100./(1+exp(-(20+160*(cR-0.5)^2)*(x_adjusted{i1,i2}(1,1:end,j1)-cR))));
% %             x_adjusted{i1,i2}(2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1) = (x_adjusted{i1,i2}(2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1)+(c/2)./(1+exp(-(20+160*(cWT-0.5)^2)*(x_adjusted{i1,i2}(2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1)-cWT))));
%        x_adjusted{i1,i2}(2,1:end,j1) = ((x_adjusted{i1,i2}(2,1:end,j1)+100*c./(1+exp(-(20+160*(cR-0.5)^2)*(x_adjusted{i1,i2}(2,1:end,j1)-cR)))) ./ (1+100*c./(1+exp(-(20+160*(cR-0.5)^2)*(1-cR)))) );
%        x_adjusted{i1,i2}(1,1:end,j1) = x_adjusted{i1,i2}(1,1:end,j1).*(1-Drug_Adherance_Uninfected{i1}(tt)); 
% %        ((x_adjusted{i1,i2}(1,1:end,j1)+(c+1e-6)./(1+exp(-(20+160*(cWT-0.5)^2)*(x_adjusted{i1,i2}(1,1:end,j1)-cWT)))) ./ (8*(c+1e-6)) );
% 
%         end
%     end
% end
  

for i1 = 1:nh
    for i2 = 1:nh
        for j1 = 1:ns
             for j2 = 1:ns
                c = c0*(1 - Drug_Adherance_Uninfected{i1}(tt)) + c1*Drug_Adherance_Uninfected{i1}(tt);
                x_adjusted{i1,i2}(1,1:end,j1) = host_x{i2}(1,1:end,j1).*(1-c)./( host_x{i2}(1,1:end,j1).*(1-c) + host_x{i2}(2,1:end,j1).*c );     
                x_adjusted{i1,i2}(2,1:end,j1) = host_x{i2}(2,1:end,j1).*c./( host_x{i2}(1,1:end,j1).*(1-c) + host_x{i2}(2,1:end,j1).*c );
             end
        end
    end
end


% a = c0*(1 - Drug_Adherance_Uninfected{i1}(tt)) + c1*Drug_Adherance_Uninfected{i1}(tt);
%                
% 
% x_adjusted{i1,i2}(1,1:end,j1) = (host_x{i2}(1,1:end,j1)-1/a*host_x{i2}(1,1:end,j1).^a)/(1-1/a);
% x_adjusted{i1,i2}(2,1:end,j1) = 1.0 - x_adjusted{i1,i2}(1,1:end,j1);

VL = cell(nh,nh);
for i = 1:nh
    for j = 1:nh
        VL{i,j} = zeros(ns,length(tt),ns);
    end
end

for i1 = 1:nh
    for i2 = 1:nh
        for j1 = 1:ns
            for j2 = 1:ns
                VL{i1,i2}(j2,1:tau1{i2}(j1),j1) = lambda1 * x_adjusted{i1,i2}(j2,1:tau1{i2}(j1),j1);
                VL{i1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1) = lambda2(j1) * x_adjusted{i1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1);
            end
        end
    end
end   


% Now assume that in the presence of the drug the infectivity of the WT has
% a maximum (due to decreased viral load). Thi is about 1% what it is for
% host not on PrEP
% for i1 = 1:nh
%     for i2 = 1:nh
%         for j1 = 1:ns
%             VL{i1,i2}(1,1:end,j1) = min( VL{i1,i2}(1,1:end,j1), 0.004*Drug_Adherance{i2}(tt) + 0.4*(1 - Drug_Adherance{i2}(tt)) );
%         end
%     end
% end


% The drug still gas big influence on WT infectivity despite the adjusted
% bit. Could also allow for drug to impact on resistant strain
  

for i1 = 1:nh
    for i2 = 1:nh
        for j1 = 1:ns
            VL{i1,i2}(1,1:end,j1) = VL{i1,i2}(1,1:end,j1) .* (1 - 0.99*Drug_Adherance_Uninfected{i1}(tt));
        end
    end
end


% Incorporate the natural mortality
for i1 = 1:nh
    for i2 = 1:nh
        for j1 = 1:ns
             for j2 = 1:ns
                VL{i1,i2}(j2,1:end,j1) = VL{i1,i2}(j2,1:end,j1) .* exp(-nu*tt);     
             end
        end
    end
end




%% Other infectivity stuff
% 
% % From the frequencies find frequency in each host, this will depend on
% % drug concentration and stuff
% 
% c0 = 0.5; % Competition parameter without the drug, when this is 1 WT
% % cannot cause infection when there is no drug
% c1 = 0.9; % Competition parameter with the drug, when this is 1 WT 
% % (susceptible) cannot cause infection in host on PrEP
% 
% x_adjusted = cell(nh,nh);
% y_adjusted = cell(nh,nh);
% for i = 1:nh
%     for j = 1:nh
%         x_adjusted{i,j} = zeros(ns,length(tt),ns);
%         y_adjusted{i,j} = zeros(ns,length(tt),ns);
%     end
% end
% % Simple way of allowing for competition, best if these don't sum to one.
% % So for recipient on the druginfectivity of each is lowered compared to
% % PrEP naive host. Transmissibility of both is affectedby the drug either
% % positively or negatively
% for i1 = 1:nh
%     for i2 = 1:nh
%         for j1 = 1:ns
%              for j2 = 1:ns
%                 c = c0*(1 - Drug_Adherance_Uninfected{i1}(tt)) + c1*Drug_Adherance_Uninfected{i1}(tt);
%                 x_adjusted{i1,i2}(1,1:end,j1) = host_x{i2}(1,1:end,j1).*(1-c)./( host_x{i2}(1,1:end,j1).*(1-c) + host_x{i2}(2,1:end,j1).*c );     
%                 x_adjusted{i1,i2}(2,1:end,j1) = host_x{i2}(2,1:end,j1).*c./( host_x{i2}(1,1:end,j1).*(1-c) + host_x{i2}(2,1:end,j1).*c );
%                 y_adjusted{i1,i2}(1,1:end,j1) = host_y{i2}(1,1:end,j1).*(1-c)./( host_y{i2}(1,1:end,j1).*(1-c) + host_y{i2}(2,1:end,j1).*c );     
%                 y_adjusted{i1,i2}(2,1:end,j1) = host_y{i2}(2,1:end,j1).*c./( host_y{i2}(1,1:end,j1).*(1-c) + host_y{i2}(2,1:end,j1).*c );
%              end
%         end
%     end
% end
% 
% VL = cell(nh,nh);
% for i = 1:nh
%     for j = 1:nh
%         VL{i,j} = zeros(ns,length(tt),ns);
%     end
% end
% 
% for i1 = 1:nh
%     for i2 = 1:nh
%         for j1 = 1:ns
%             for j2 = 1:ns
%                 VL{i1,i2}(j2,1:tau1{i2}(j1),j1) = lambda1 * ( (1-transmission) * x_adjusted{i1,i2}(j2,1:tau1{i2}(j1),j1) + transmission * y_adjusted{i1,i2}(j2,1:tau1{i2}(j1),j1));
%                 VL{i1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1) = lambda2(j1) * ( (1-transmission) * x_adjusted{i1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1) + transmission * y_adjusted{i1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1));
%             end
%         end
%     end
% end   
% 
% 
% % Now assume that in the presence of the drug the infectivity of the WT has
% % a maximum (due to decreased viral load). Thi is about 1% what it is for
% % host not on PrEP
% for i1 = 1:nh
%     for i2 = 1:nh
%         for j1 = 1:ns
%             VL{i1,i2}(1,1:end,j1) = min( VL{i1,i2}(1,1:end,j1), 0.004*Drug_Adherance{i2}(tt) + 0.4*(1 - Drug_Adherance{i2}(tt)) );
%         end
%     end
% end
% 
% 
% % The drug still gas big influence on WT infectivity despite the adjusted
% % bit. Could also allow for drug to impact on resistant strain
%   
% 
% for i1 = 1:nh
%     for i2 = 1:nh
%         for j1 = 1:ns
%             VL{i1,i2}(1,1:end,j1) = VL{i1,i2}(1,1:end,j1) .* (1 - 0.99*Drug_Adherance_Uninfected{i1}(tt));
%         end
%     end
% end
% 
% 
% % Incorporate the natural mortality
% for i1 = 1:nh
%     for i2 = 1:nh
%         for j1 = 1:ns
%              for j2 = 1:ns
%                 VL{i1,i2}(j2,1:end,j1) = VL{i1,i2}(j2,1:end,j1) .* exp(-nu*tt);     
%              end
%         end
%     end
% end
% 
% % % The viral load of the WT for someone on PrEP will be very low, but for
% % % resistant strain it will be about normal?(given the assumptions about 
% % % fitness, i.e. that the drug does not affect the resistant strain).
% % 
% % for i1 = 1:nh
% %     for i2 = 1:nh
% %         for j1 = 1:ns
% %             VL{i1,i2}(1,1:end,j1) = VL{i1,i2}(1,1:end,j1) .* (1 - Drug_Adherance{i2}(tt));
% %         end
% %     end
% % end  
% % 
% % % This is now the infectivity, to account for competition between strains 
% % % during transmission the infectivity (or maybe frequency is adjusted?)
% % VL_adjusted = cell(nh,nh);
% % for i = 1:nh
% %     for j = 1:nh
% %         VL_adjusted{i,j} = zeros(ns,length(tt),ns);
% %     end
% % end
% % for i1 = 1:nh
% %     for i2 = 1:nh
% %         for j1 = 1:ns
% %              for j3 = 1:ns
% %                 VL_adjusted{i1,i2}(j3,1:end,j1) = VL{i1,i2}(j3,1:end,j1);     
% %              end
% %         end
% %     end
% % end
% % 
% % cs = 0.5;
% % for i1 = 1:nh
% %     for j1 = 1:ns
% %         VL_adjusted{2,i1}(1,1:end,j1) = VL{2,i1}(1,1:end,j1)*(1-cs) ./ ( VL{2,i1}(1,1:end,j1)*(1-cs) + VL{2,i1}(2,1:end,j1)*cs );
% %         VL_adjusted{2,i1}(2,1:end,j1) = VL{2,i1}(2,1:end,j1)*cs ./ ( VL{2,i1}(1,1:end,j1)*(1-cs) + VL{2,i1}(2,1:end,j1)*cs );
% %     end
% % end
% % 
% % for i1 = 1:nh
% %     for i2 = 1:nh
% %         for j1 = 1:ns
% %             VL{i1,i2}(1,1:end,j1) =  VL{i1,i2}(1,1:end,j1) .* (1 - Drug_Adherance_Uninfected{i1}(tt));
% %         end
% %     end
% % end
% 
% 
% 
% 
% %%%%%%
% %%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% %%%%%%%%% The simpler version with VL x 1 or VL x 0, for PrEP and not PrEP
% %%%%%%%%%%%%%%%%%%%%
% % 
% % VL = cell(nh,nh);
% % for i = 1:nh
% %     for j = 1:nh
% %         VL{i,j} = zeros(ns,length(tt),ns);
% %     end
% % end
% % 
% % for i1 = 1:nh
% %     for i2 = 1:nh
% %         for j1 = 1:ns
% %             for j2 = 1:ns
% %                 VL{i1,i2}(j2,1:tau0{i2}(j1),j1) = lambda0;
% %                 VL{i1,i2}(j2,(tau0{i2}(j1)+1):tau1{i2}(j1),j1) = lambda1;
% %                 VL{i1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1) = lambda2(j1);
% %             end
% %         end
% %     end
% % end   
% % 
% % % not on PrEP chance of being infected donor frequency times VL
% % for i2 = 1:nh
% %     for j1 = 1:ns
% %         for j2 = 1:ns
% %             VL{1,i2}(j2,1:tau0{i2}(j1),j1) = VL{1,i2}(j2,1:tau0{i2}(j1),j1) .* host_x{i2}(j2,1:tau0{i2}(j1),j1);
% %             VL{1,i2}(j2,(tau0{i2}(j1)+1):tau1{i2}(j1),j1) = VL{1,i2}(j2,(tau0{i2}(j1)+1):tau1{i2}(j1),j1) .* host_x{i2}(j2,(tau0{i2}(j1)+1):tau1{i2}(j1),j1);
% %             VL{1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1) = VL{1,i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1) .* host_x{i2}(j2,(tau1{i2}(j1)+1):tau2{i2}(j1),j1);
% %         end
% %     end
% % end
% % 
% % % ON PrEP chance of being infected with resistant strain 1 and
% % % non-resistant 0
% % 
% % for i1 = 2:nh
% %     for i2 = 1:nh
% %         for j1 = 1:ns
% %             VL{i1,i2}(1,1:end,j1) =  VL{i1,i2}(1,1:end,j1) .* (1 - Drug_Adherance_Uninfected{i1}(tt));
% %             VL{i1,i2}(2,1:end,j1) =  VL{i1,i2}(2,1:end,j1) .* Drug_Adherance_Uninfected{i1}(tt); 
% %         end
% %     end
% % end
% % 
% % 
% % 
% % % Incorporate the natural mortality
% % for i1 = 1:nh
% %     for i2 = 1:nh
% %         for j1 = 1:ns
% %              for j3 = 1:ns
% %                 VL{i1,i2}(j3,1:end,j1) = VL{i1,i2}(j3,1:end,j1) .* exp(-nu*tt);     
% %              end
% %         end
% %     end
% % end
% % 
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %     
% %     
% %     
% %     
% %     
% % % Strain 1 VL is reduced by the drug, so it is ~0 when a lot of drug, when
% % % initiated by either strain. It won't go completely to 0 (andif it did
% % % there would be NaNs in the solution). When infected initiLLY WITH STRAIN
% % % 2, the VL follows the dynnamics in the absence of drug
% % % for i = 1:nh
% % %     for j = 1:nh
% % %        VL{i,j}(1,1:end,1) =  max(VL{i,j}(1,1:end,1) .* (1 - Drug_Adherance{j}(tt)), 0.001);
% % %     end
% % % end
% % % % Lets say after 30 days the VL of the resistant strain reaches the spvl,
% % % % since there is definitely some present 
% % % tStart = find(Drug_Adherance{i}(tt) >= 0, 1);
% % % tStop = find(tt >= tt(tStart) + 0.0822, 1);
% % % 
% % % for i = 1:nh
% % %     for j = 1:nh
% % %         VL{i,j}(2,tStart:tStop,1) =  min(VL{i,j}(2,tStart:tStop,1) .* (1 - Drug_Adherance{j}(tt(tStart:tStop))), 0.001);
% % %      end
% % % end
% % 
% % % Make transmission from host 1 to host 2 almost zero, this slow dynamics
% % % maybe? Should be quite low anyway cause of the low freq
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % These are essentially the x(1-P) bits, but it is assumed that P = 1 for
% % % all time unless the person gets infected
% %     VL{2,1}(2,1:end,1) =  min(VL{2,1}(2,1:end,1), 1e-20);
% %     VL{2,1}(2,1:end,2) =  min(VL{2,1}(2,1:end,2), 1e-20);
% % % Host 1 (not on PrEP) should not be able to transmit virus to someone on
% % % PrEP (because host 2 is on PrEP)
% %     VL{2,1}(1,1:end,1) =  min(VL{2,1}(1,1:end,1), 1e-20);
% %     VL{2,1}(1,1:end,2) =  min(VL{2,1}(1,1:end,2), 1e-20);
% % % If host 2 stops PrEP after getting the virus then WT increases but we 
% % % still assume that the hosts who don't have the virus remain on PrEP
% %     VL{2,2}(1,1:end,1) =  min(VL{2,2}(1,1:end,1), 1e-20);
% %     VL{2,2}(1,1:end,2) =  min(VL{2,2}(1,1:end,2), 1e-20);
% %  
% %     
% %     % How infectious WT is to host 3 depends on how much drug they have in
% %     % system 
% %     VL{3,1}(1,1:end,1) = max(VL{3,1}(1,1:end,1) .* (1 - Drug_Adherance{3}(tt)), 0.001);
% %     VL{3,1}(1,1:end,2) = max(VL{3,1}(1,1:end,2) .* (1 - Drug_Adherance{3}(tt)), 0.001);
% %  
% %     VL{3,2}(1,1:end,1) = max(VL{3,2}(1,1:end,1) .* (1 - Drug_Adherance{3}(tt)), 0.001);
% %     VL{3,2}(1,1:end,2) = max(VL{3,2}(1,1:end,2) .* (1 - Drug_Adherance{3}(tt)), 0.001);
% %     
% %     
% % % Does this need to be chnged so that resistant is less infectious when they off the drug 
% % % VL{3,2}(:,1:length(tt),2)
% % 
% % % We assume the infectvity of the WT strain (1) is further altered by the
% % % concentration of drug in the recipient host, using the same (1-P) thing.
% % % So that when recipient has maximum anount of  PrEP they cannot be
% % % infected by the WT.The max is so that there are no zeros in the next
% % % generation matrix
% % 
% % 
% % for i = 1:nh
% %     for j = 1:nh
% %        VL{i,j}(1,1:end,1) =  max(VL{i,j}(1,1:end,1) .* (1 - Drug_Adherance{i}(tt)), 0.001);
% %        VL{i,j}(1,1:end,2) =  max(VL{i,j}(1,1:end,2) .* (1 - Drug_Adherance{i}(tt)), 0.001);
% %     end
% % end
% % 
% % 
% % 
% % % Incorporate the natural mortality
% % for i1 = 1:nh
% %     for i2 = 1:nh
% %         for j1 = 1:ns
% %              for j3 = 1:ns
% %                 VL{i1,i2}(j3,1:end,j1) = VL{i1,i2}(j3,1:end,j1) .* exp(-nu*tt);     
% %              end
% %         end
% %     end
% % end
% 
% 
% %%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%

%%
%%% BETWEEN-HOST DYNAMICS %%%

%%% Initialisation %%%

ICint = appTmax; % Duration of interval of initial conditions, before tmin, during which incidence is constant at I0
iT0 = ICint/dt; % Incidence is 0 up to t(iT0)

t = (tmin-ICint):dt:tmax;
lt = length(t);

host_incidence = cell(nh,nh);

setup_1 = cell(nh,1);
setup_ns = cell(nh,1);
for i = 1:nh
    setup_1{i} = zeros(1,lt);
    setup_ns{i} = zeros(ns,lt);
    for j = 1:nh
        host_incidence{i,j} = zeros(ns,lt);
    end
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
I0_2 = cell(nh,nh);
for i = 1:nh
    I0{i} = zeros(ns,1);   % Column vector of initial values for each strain
    for j = 1:nh
        I0_2{i,j} = zeros(ns,1); % Initial condition for looking at hosts individually
    end
end

I0{starthost}(startstrain) = 1;          % Initial condition for the incidence in the population, with only strain 1 present, we assume it starts in someone not on PrEP
I0_2{starthost,1}(startstrain) = 1;        % Assume host 1 caused initial infection
S0 = N0;

for iic = 1:iT0
    for st = 1:nh
        I{st}(:,iic) = I0{st};
        Isum = 0;%zeros(iT0,1);
        for i = 1:nh 
            Isum = Isum + sum(I0{i});
            host_incidence{st,i}(:,iic) = I0_2{st,i};
        end
        Ibar{st}(:,iic) = I0{st} / Isum; % maybe don't need totals in different hosts
        Itot{st}(iic) = sum(I0{st});
        J{st}(:,iic) = zeros(ns,1);         % Could do better, probably
        Jbar{st}(:,iic) = zeros(ns,1); 
        Jtot{st}(iic) = 0; 
        S{st}(iic) = Hprobs(st)*S0; 
    end
    N(iic) = N0;
end

N(iT0) = N0;
for i = 1:nh
    S{i}(iT0) = Hprobs(i) * (N(iT0) - Isum * dt * trapz( disc_vec ));
end
%     
% S{1}(iT0) = Hprops(1) * (N(iT0) - (Itot{1}(iT0)+Itot{2}(iT0)) * dt * trapz( disc_vec ));
% S{2}(iT0) = Hprops(2) * (N(iT0) - (Itot{1}(iT0)+Itot{2}(iT0)) * dt * trapz( disc_vec ));


% Solve between-host system
FOI = cell(nh,1);         % Force of infection of different strains over time for different hosts
FOI_temp2 = cell(nh,1);
FOI_host2 = cell(nh,nh);
for i = 1:nh
    FOI{i} = zeros(ns,lt);
    FOI_temp2{i} = zeros(ns,1);
    for j = 1:nh
       FOI_host2{i,j} = zeros(ns,1);
    end
end


for it = (iT0+1):lt
    discI = cell(nh,1);
    for discI_ind = 1:nh
        discI{discI_ind} = zeros(ns,1);
    end
    for i = 1:ns % The strain that initiates infection in new host 
        FOI_temp = cell(nh,1);  % FOI from type j at any point in the past (acting on i, but not depending explicitly, i.e. overwritten every loop)
        FOI_host = cell(nh,nh); % To see what host are causing what infections
        for FOI_temp_ind = 1:nh
            FOI_temp{FOI_temp_ind} = zeros(ns,ltt);
            for FOI_host_ind = 1:nh
                FOI_host{FOI_temp_ind,FOI_host_ind} = zeros(ns,ltt);
            end
        end       
        for j = 1:ns % The strain that initiated infection in the donor (determines VL)
            for FOI_temp_ind = 1:nh
                FOI_temp{FOI_temp_ind}(j,1) = 0;    % The first point of this vector is always 0, because it refers to the time point t(it), the one we are calculating the solution at
            end
            for FOI_temp_ind1 = 1:nh
               for FOI_temp_ind2 = 1:nh
                    for itt = 2:min(it,tau2{FOI_temp_ind2}(j))
                   % The force of infection of virus initiated by strain j for each host
                        FOI_temp{FOI_temp_ind1}(j,itt) = FOI_temp{FOI_temp_ind1}(j,itt) + VL{FOI_temp_ind1,FOI_temp_ind2}(i,itt,j) * I{FOI_temp_ind2}(j,it-itt+1);% host_x{FOI_temp_ind2}(i,itt,j) * I{FOI_temp_ind2}(j,it-itt+1);
                   % The force of infection due to each host and strain,
                   % i.e. each term in the sum
                        FOI_host{FOI_temp_ind1,FOI_temp_ind2}(j,itt) = VL{FOI_temp_ind1,FOI_temp_ind2}(i,itt,j) * I{FOI_temp_ind2}(j,it-itt+1);% host_x{FOI_temp_ind2}(i,itt,j) * I{FOI_temp_ind2}(j,it-itt+1);
                    end
                    
               end
            end
                
            for FOI_temp2_ind = 1:nh
                FOI_temp2{FOI_temp2_ind}(j) = dt * trapz( FOI_temp{FOI_temp2_ind}(j,1:min(it,tau2{FOI_temp2_ind}(j))) );     % Contibution of j to the FOI acting on i
                for FOI_host2_ind = 1:nh
                    FOI_host2{FOI_temp2_ind,FOI_host2_ind}(j) = dt * trapz( FOI_host{FOI_temp2_ind,FOI_host2_ind}(j,1:min(it,tau2{FOI_host2_ind}(j))) );     % Contibution of j to the FOI acting on i
                end
            end
        end
        for host = 1:nh
            for ihost = 1:nh
                host_incidence{host,ihost}(i,it) = S{host}(it-1) / N(it-1) * sum(FOI_host2{host,ihost}); % Sum over strains   
            end
            FOI{host}(i,it) = sum( FOI_temp2{host} );%   FOI_temp2{host}(i);%      should this be a sum? It was already summed before the integral
            I{host}(i,it) = S{host}(it-1) * FOI{host}(i,it) / N(it-1);
            J{host}(i,it) = dt * trapz( I{host}(i,max(1,it:-1:it-tau2{host}(i)+1)) .* disc_vec(1:tau2{host}(i)) );
            discI{host}(i) = I{host}(i,max(1,it-tau2{host}(i))) * disc_last{host}(i);
        end      
% discI is the sum of incidence bit in the equation for N, this goes in the
% Euler step bit. So for N(t) we need this at t-dt, so doesthis mean the
% max works?
        
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
        Hprobs = Hprobs2;
    end
    
    N(it) = N(it-1) + dt * ( B - nu * N(it-1) - discIsum );
    
    Sus = N(it) - Jsum;
    for host = 1:nh
%%%%%%%%% For fixed proportion use this
        S{host}(it) = Hprobs(host) * Sus;
%%%%%%%%%
%%%%%%%%% For a fixed "birth" rate use this
% Euler step perhaps not the best way to do this? Can't use Runge
% Kutta beacuse we would need to evaluate the RHS at dt/2, this is not
% calculated for incidence etc.
% Need derivative of infected 
%         Idiff = gradient(sum(J{host}));
%         S{host}(it) = S{host}(it-1) + dt * ( Hprobs(host) * B - nu * ( S{host}(it-1) + 0*sum(J{host}(:,it-1)) ) - 0*sum(discI{host}) - Idiff(it-1) );
%%%%%%%%% 
        Ibar{host}(:,it) = I{host}(:,it) / Isum;
        Jbar{host}(:,it) = J{host}(:,it) / Jsum;
    end

end

%% FIGURES 

Prevalence = Jbar{1};
Susceptibles = S{1};
Infected = sum(J{1});
for i = 2:nh
    Prevalence = Prevalence + Jbar{i};
    Susceptibles = Susceptibles + S{i};
    Infected = Infected + sum(J{i});
end




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
ylabel('Prevalence','Fontsize',20,'interpreter','latex');
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







%% The R_0, the different infectivities can be put into a single matrix

for i = 0:nh*ns-1
    for j = 0:nh*ns-1
        NGM(i+1,j+1) = Hprobs( mod(floor(i/ns),nh) + 1 ) * dt * trapz( VL{mod(floor(i/ns),nh) + 1, mod(floor(j/ns),nh) + 1}((mod(i,ns) + 1),:,(mod(j,ns) + 1)) );
   end
end

 time_vector = [];

for i = 1:nh
    for j = 1:ns
        time_vector = [time_vector; tt(tau2{i}(j))];
    end
end



     [V,D] = eig(NGM);
     
    [maxD,ind] = max(D(:));
    [m,n] = ind2sub(size(D),ind);

     R_0 = maxD;
     % Should be the same as (or close to)
     S_total = 0;
     for i = 1:nh
         S_total = S_total + S{i}(end);
     end
     R_0a =  N(end)/(S_total);
     
     
     eql = V(:,n)/sum(V(:,n));
     % Incidence at equilibrium
     Ifac = (B * (R_0 - 1))/(R_0 - sum(eql.*exp(-nu*time_vector)));
     Ieql = Ifac * eql;
     
     Ne = ( B - sum(Ieql.*exp(-nu*time_vector)) )/nu;
     
     IeqlNum = [];
     for i = 1:nh
         for j = 1:ns
             IeqlNum = [IeqlNum; I{i}(j,end)];
         end
     end
     
     eqlNum = IeqlNum/sum(IeqlNum);
     
     ResPrev = sum(eqlNum(2:2:end));
     
     infectedNum = [];
     for i = 1:nh
         for j = 1:ns
             infectedNum = [infectedNum; J{i}(j,end)];
         end
     end
     
     infectedAn = Ieql.*(1-exp(-nu*time_vector))/nu;
   




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
fr = 1;
figure
plot(t, (host_incidence{1,fr}(:,1:end)+host_incidence{2,fr}(:,1:end)+host_incidence{3,fr}(:,1:end)+host_incidence{4,fr}(:,1:end)) / Hprobs(fr))
 
 