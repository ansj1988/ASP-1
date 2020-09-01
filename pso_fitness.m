function fitness = pso_fitness(population,signal,noise,test)

% Recovering the swarm size
swarm_size = length(population(:,1));

for k = 1:swarm_size

    % Recovering the particles coordinates
    stepsize = population(k,1);
    order = population(k,2);
    
    if test == 0
        
        % Evaluating the performance of each particle
        e = nlms_evaluation(signal,noise,order,stepsize);
        
        % Classical MSE fitness
%         fitness(k) = mean(e(:,1).^2);                                                     
        
        % Proposed fitness
        vec1 = e(end-4999:end,1);
        vec2 = noise(end-4999:end);     
        fitness(k) = abs(corr(vec1,vec2','Type','Pearson'));
        
    % Ackley function
    elseif test == 1
        fitness(k) = -20*exp(-0.2*sqrt(0.5*(stepsize^2+order^2))) - exp(sqrt(0.5*(cos(2*pi*stepsize)+cos(2*pi*order)))) + 20 + exp(1); 
    end
    
end
end