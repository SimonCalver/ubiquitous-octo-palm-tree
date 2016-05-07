function [x,y] = solve_qsODE_PrEP(tau,z0,ns,g,fmut,bmut,rmut,k,a,rL)
% This function solves the within-host quasispecies equation, using
% matlab's build in ODE solver ode45. The initial vector z0 consists of the
% initial strain distribution in the active compartment (x0), followed by
% the initial strain distribution in the reservoir (y0). I.e. z0 =
% (x0,y0). Other parameters:
% ns    number of strains
% g     vector with strain-specific replication rates
% fmut  probability of a forward mutation (strain i -> i+1)
% bmut  probability of a backward mutation (strain i -> i-1)
% rmut  probability of a random other mutation
% k     probability that newly infected cells become latent
% a     activation rate of latently infected cells
% rL    relative size of the reservoir

% Check if number of strains corresponds with length of initial vector
if length(z0) ~= 2*ns
    disp('Error solve_qsODE1: initial vector should have length 2ns')
    return
end

% Solve the full ODE (including chemostat)
[t,z] = ode45(@RHS,tau,z0);

%options = odeset('RelTol',1e-2);
%[t,z] = ode45(@RHS,tau,z0,options);

x = (z(:,1:ns))';
y = (z(:,(ns+1):2*ns))';

    % Construct the RHS of the ODE-system
    function fz = RHS(t,z)
        fz = zeros(2*ns,1);
        
        % Construction of full matrix Q, including reproduction/mutation
        % and inflow/outflow to/from reservoir
        Q = zeros(2*ns);

        %Mutation matrix
        M = rmut*ones(ns,ns);
        M(2,1) = fmut;                 % Set border values in matrix
        M((ns-1),ns) = bmut;
        M(1,1) = (1-fmut-(ns-2)*rmut);
        M(ns,ns) = (1-bmut-(ns-2)*rmut);
        for j = 2:(ns-1)              % Set other values (direct mutations + prob of no mutation on diagonal)
            M((j-1),j) = bmut;
            M(j,j) = (1-fmut-bmut-(ns-3)*rmut);
            M((j+1),j) = fmut;
        end

        % Submatrices
        R = M*diag(g(t));              % Reproduction-mutation
        A = a*eye(ns);              % Leaving reservoir

        % Composition of full matrix Q from submatrices
        Q(1:ns, 1:ns) = (1-k)*R;
        Q((ns+1):2*ns, 1:ns) = (k/rL)*R;
        Q(1:ns, (ns+1):2*ns) = rL*A;
        Q((ns+1):2*ns, (ns+1):2*ns) = -A;

        % Construction of chemostat matrix
        totgrowth = (z(1:ns))'*g(t);
        diagvec = zeros(2*ns,1);
        diagvec(1:ns) = ((1-k)*totgrowth + rL*a)*ones(ns,1);
        diagvec((ns+1):2*ns) = ((k/rL)*totgrowth-a)*ones(ns,1);
        C = diag(diagvec);
        
        fz = Q*z - C*z;             % Return value of the RHS of the ode (the value of f(t,z))
    end

end