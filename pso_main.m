function [optimum_stepsize,optimum_order] = pso_main(signal,noise)

%%% COMMENTS %%%
% 1. Classical PSO is for maximization; this is a minimization version
% 2. If you want to check the consistency of this code, please call the 'pso_main' function 
% using any inputs and change the 'test' variable to 1. 
% In this case, a test function - Ackley function - with minimum located in
% (0,0) is used.
test = 0;

%% Step 1: setting the main parameters

% Boundaries of the search space
if test == 0  
    stepsize_min = 0.01;    % Minimum stepsize
    stepsize_max = 0.1;     % Maximum stepsize
    order_min = 10;         % Minimum order
    order_max = 50;         % Maximum order
elseif test == 1    
    stepsize_min = -40;     % Minimum stepsize
    stepsize_max = 40;      % Maximum stepsize
    order_min = -40;        % Minimum order
    order_max = 40;         % Maximum order
end

% Population
M = 50;                 % Swarm size
tam_sub = 5;            % Sub-swarm size
subswarm = 1;           % Enables sub-swarms
N_sub = M/tam_sub;      % Number of sub-swarms

if test==1
    subswarm=0;
end

% Dynamics of the process
acc1 = 2;               % Acceleration constant
acc2 = 2;               % Acceleration constant
vel = zeros(M,2);       % Velocities

% Convergence criteria
fitness_min = -1e-8;     % Lower limit of the fitness function
iterations_max = 100;    % Maximum number of iterations

%% Step 2: first swarm

% Generating the initial population
pso_iteration = 1;                                                                 % Counter
population = pso_first_swarm(M,order_min,order_max,stepsize_min,stepsize_max); % Coordinates of all particles 

% Graph
% First generation of particles over the search space
figure(1)
for k = 1:M
    plot(population(k,1),population(k,2),'*','color','r')
    xlabel('Step-size (\mu)')
    ylabel('Order (L)')
    xlim([stepsize_min stepsize_max])
    ylim([order_min order_max])
    hold on
end

% Evaluating the fitness values
fitness = pso_fitness(population,signal,noise,test);

% Storing the coordinates and fitness values of each particle
memory_fitness(:,pso_iteration) = fitness;
memory_stepsize(:,pso_iteration) = population(:,1);
memory_order(:,pso_iteration) = population(:,2);
    
% Calculating gbest
[gbest(3),gbest_id] = min(fitness);               % Particle with the minimum fitness of the swarm
gbest(1) = memory_stepsize(gbest_id,pso_iteration);   % gbest x-coordinate (step-size)
gbest(2) = memory_order(gbest_id,pso_iteration);      % gbest y-coordinate (order)

% Storing the gbest coordinates and fitness value
memory_gbest(1,pso_iteration) = gbest(1); % gbest x-coordinate
memory_gbest(2,pso_iteration) = gbest(2); % gbest y-coordinate
memory_gbest(3,pso_iteration) = gbest(3); % gbest fitness
 
% Calculating the local gbest (only if you are using sub-swarms)
if subswarm == 1
    for k=0:N_sub-1
        [gbest_local(3,k+1),gbestlocal_id] = min(fitness(k*tam_sub+1:(k+1)*tam_sub));   % Particle with the minimum fitness of the sub-swarm
        gbestlocal_id = k*tam_sub + gbestlocal_id;                                      % Converting (1,N_sub) scale to (1,M) range
        gbest_local(1,k+1) = memory_stepsize(gbestlocal_id,pso_iteration);                  % Local gbest x-coordinate
        gbest_local(2,k+1) = memory_order(gbestlocal_id,pso_iteration);                     % Local gbest y-coordinate
    end
else
    gbest_local = 0;
end
    
% Calculating pbest 
% In the first swarm, there is no historical data from past iterations so pbest equals the 
% particles initial coordinates and fitness
for k=1:M
    pbest(k,1) = population(k,1); % pbest x-coordinate
    pbest(k,2) = population(k,2); % pbest y-coordinate
    pbest(k,3) = fitness(k);      % pbest fitness
end

%% Step 3: updating the swarms until convergence

% Decision variable:
% decision = 1: process continues
% decision = 0: process stops
decision = 1; 

while (decision == 1)

    pso_iteration = pso_iteration + 1
    [population, vel] = pso_new_swarm(pso_iteration,population,gbest,pbest,acc1,acc2,vel,M,order_min,order_max,stepsize_min,stepsize_max,iterations_max,gbest_local,tam_sub,subswarm);
    fitness = pso_fitness(population,signal,noise,test);
    
    % If you want to display the particles distribution over the search space 
    % at every iteration, please enable the following lines
    
    % figure
    % for k = 1:M
    % plot(population(k,1),population(k,2),'*','color','k')
    % xlabel('Step-size (\mu)')
    % ylabel('Order (L)')
    % xlim([stepsize_min stepsize_max])
    % ylim([order_min order_max])
    % hold on
    % end
    
    % Storing the coordinates and fitness values of each particle
    memory_fitness(:,pso_iteration) = fitness;
    memory_stepsize(:,pso_iteration) = population(:,1);
    memory_order(:,pso_iteration) = population(:,2);

    % Calculating gbest
    [gbest(3),gbest_id] = min(fitness);               % Particle with the minimum fitness of the swarm
    gbest(1) = memory_stepsize(gbest_id,pso_iteration);   % gbest x-coordinate (step-size)
    gbest(2) = memory_order(gbest_id,pso_iteration);      % gbest y-coordinate (order)

    % Storing the gbest coordinates and fitness value
    memory_gbest(1,pso_iteration) = gbest(1); % gbest x-coordinate
    memory_gbest(2,pso_iteration) = gbest(2); % gbest y-coordinate
    memory_gbest(3,pso_iteration) = gbest(3); % gbest fitness
    
    % Calculating the local gbest (only if you are using sub-swarms)
    if subswarm == 1
        for k=0:N_sub-1
            [gbest_local(1,3,k+1),gbestlocal_id] = min(fitness(k*tam_sub+1:(k+1)*tam_sub)); % Particle with the minimum fitness of the sub-swarm
            gbestlocal_id = k*tam_sub + gbestlocal_id;                                      % Converting (1,N_sub) scale to (1,M) range
            gbest_local(1,1,k+1) = memory_stepsize(gbestlocal_id,pso_iteration);                % Local gbest x-coordinate
            gbest_local(1,2,k+1) = memory_order(gbestlocal_id,pso_iteration);                   % Local gbest y-coordinate
        end
    else
        gbest_local = 0;
    end

    % Calculating pbest 
    for k=1:M     
        pbest(k,1) = population(k,1);          % pbest x-coordinate
        pbest(k,2) = population(k,2);          % pbest y-coordinate
        pbest(k,3) = min(memory_fitness(k,:)); % Minimum fitness of the kth particle along the iterations   
    end

    % Evaluating if the process stops or continues
    if pso_iteration >= iterations_max  % 1st criterium: maximum number of iterations
        decision = 0; 
    else                            % 2nd criterium: lower limit of the fitness function
        for k = 1:M      
            if abs(fitness(k)) < abs(fitness_min)  
                fitness(k)
                decision = 0;
            end
        end
    end
       
end

% The optimum parameters are defined as the gbest coordinates from the last
% iteration
optimum_stepsize = memory_gbest(1,end);
optimum_order = memory_gbest(2,end);

%% Graphs

% Evolution of the gbest fitness along the iterations
figure
plot(abs(memory_gbest(3,:)))
ylabel('Fitness')
xlabel('Iterations')
xlim([0 iterations_max])

% Last generation of particles over the search space
figure(1)
hold on
for k = 1:M
    plot(population(k,1),population(k,2),'*','color','b')
    xlim([stepsize_min stepsize_max])
    ylim([order_min order_max])
    hold on
end

end

