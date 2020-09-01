function population = pso_first_swarm(M,order_min,order_max,stepsize_min,stepsize_max)

% This variable will store the coordinates of the swarm particles
population = zeros(M,2); 

for k = 1:M
    
    stepsize = stepsize_min + (stepsize_max-stepsize_min)*rand; % Uniform distribution from stepsize_min to stepsize_max
    order = round(order_min + (order_max-order_min)*rand);      % Uniform distribution from order_min to order_max (integer numbers only) 
    population(k,1) = stepsize;                                 % Stores the stepsize of each particle
    population(k,2) = order;                                    % Stores the order of each particle
                                
end

end