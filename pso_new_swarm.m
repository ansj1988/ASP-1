function [population, vel] = pso_new_swarm(iteration,population,gbest,pbest,acc1,acc2,vel,M,order_min,order_max,stepsize_min,stepsize_max,iterations_max,gbest_local,tam_sub,subswarm)

% Some examples of inertial weight functions
% inertial_weigth = 0.9-(0.4/iterations_max)*iteration;   % Initial velocities: low, final vel.: high
% inertial_weigth = 0.9;                                  % Velocities don't change significantly 
inertial_weigth = 0.5+(0.4/iterations_max)*iteration;   % Initial velocities: high, final vel.: low

% Velocity limits
vmax(1) = (0.5/2)*(stepsize_max-stepsize_min);
vmax(2) = (0.5/2)*(order_max-order_min);

% Counter used to identify the sub-swarms
cont = 0;

% Updating velocities and positions
for k=1:M
    
    % With sub-swarms
    if subswarm == 1
    
        vel(k,1) = inertial_weigth*vel(k,1) + acc1*rand*(gbest_local(1,cont+1)-population(k,1))' + acc2*rand*(pbest(k,1)-population(k,1))';
        vel(k,2) = inertial_weigth*vel(k,2) + acc1*rand*(gbest_local(2,cont+1)-population(k,2))' + acc2*rand*(pbest(k,2)-population(k,2))';

        % Changing from a sub-swarm to another
        if rem(k,tam_sub) == 0
            cont = cont + 1;
        end

    % No sub-swarms
    elseif subswarm == 0
   
        vel(k,1) = inertial_weigth*vel(k,1) + acc1*rand*(gbest(1)-population(k,1))' + acc2*rand*(pbest(k,1)-population(k,1))';
        vel(k,2) = inertial_weigth*vel(k,2) + acc1*rand*(gbest(2)-population(k,2))' + acc2*rand*(pbest(k,2)-population(k,2))';  
   
    end

    % Checking the feasibility of the velocities
    if vel(k,1) > vmax(1)
        vel(k,1) = vmax(1);
    elseif vel(k,1) < -vmax(1)
        vel(k,1) = -vmax(1);
    end
    
    % Checking the feasibility of the velocities
    if vel(k,2) > vmax(2)
        vel(k,2) = vmax(2);
    elseif vel(k,2) < -vmax(2)
        vel(k,2) = -vmax(2);
    end
    
    % Updating the positions of the particles on the search space
    population(k,1) = population(k,1) + vel(k,1);        % x-coordinate
    population(k,2) = round(population(k,2) + vel(k,2)); % y-coordinate
    
    % Checking the feasibility of the positions
    if population(k,1) < stepsize_min
        population(k,1) = stepsize_min;
    end
    if population(k,2) < order_min
        population(k,2) = order_min;
    end
    
    % Checking the feasibility of the positions
    if population(k,1) > stepsize_max
        population(k,1) = stepsize_max;
    end
    if population(k,2) > order_max
        population(k,2) = order_max;
    end
    
end

end