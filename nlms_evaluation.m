function e = nlms_evaluation(signal,noise,order,stepsize)

% Initializing local variables
size = length(signal);     % Signal size
d = zeros(1,size);         % Regressor vector
u = zeros(1,order);        % Input vector
w = zeros(order,size);     % 1st dimension: weights indexes, 2nd dim: time
e = zeros(size,1);         % Error vectors
epsilon = 0.00001;         % Regularization parameter

% Calculations
for i=order:size                                                                 % Time index
    u = wrev(noise(i-order+1:i));                                                % Builds the M-size regressor vector
    d(i) = signal(i);                                                            % Desired signal            
    e(i) = d(i)-u*w(:,i-1);                                                      % Error signal
    nlms_stepsize = stepsize/(epsilon + var(u)*order);                           % Normalized step-size
    w(1:order,i) = w(1:order,i-1) + nlms_stepsize*ctranspose(u(1:order))*e(i,1); % Updating the weights
end

end