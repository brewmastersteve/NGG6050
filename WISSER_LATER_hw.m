% ChatGPT conversation used to produce this code, some modifications were
% made based on error messages

% https://chat.openai.com/share/4a8651fe-f15d-4e7b-b89d-75cb47015afd
% 
% Get the data
%   
%   Use this code to get a data set (array of RTs from a single condition) 
%   to fit, already preprocessed to include correct trials only and remove
%   outliers (including express saccades). See later_getData for details
data = later_getData([], [], 0.2);
RTs = data{1};
clear data

%lets first create the 1/RT values to "gaussianize" the data

% Create a variable "adjusted_RTs" of the same size as "RTs"
adjusted_RTs = zeros(size(RTs));

% Loop through each element of "RTs" and apply the formula
for i = 1:length(RTs)
    adjusted_RTs(i) = 1 / RTs(i);
end

% Calculate the mean and standard deviation of the "gaussian" data
meanRT = mean(adjusted_RTs);
stdDevRT = std(adjusted_RTs);

% Create a parameter array "params"
params = [meanRT, stdDevRT];


%Now we're ready to define the objective function
%using a function defined below calculate negative sum log likelihoods
laterErrFcn = @(fits) -sum(logLikelihood(fits, params(1), params(2)));



%define initial bounds
lowerBounds = [0.001, 0.001];
upperBounds = [1000, 1000];

% Choose initial conditions within the bounds
initial_mu_guess = meanRT;    
initial_sigma_guess = stdDevRT; 

% Create an array for initial conditions
initialValues = [initial_mu_guess, initial_sigma_guess];



%Run fits
opts = optimoptions(@fmincon,    ... % "function minimization with constraints"
   'Algorithm',   'active-set',  ...
   'MaxIter',     3000,          ...
   'MaxFunEvals', 3000);

% Definine the "optimization problem" using variables defined above
problem = createOptimProblem('fmincon',    ...
    'objective',   laterErrFcn,     ... % Use the objective function
    'x0',          initialValues,   ... % Initial conditions
    'lb',          lowerBounds,     ... % Parameter lower bounds
    'ub',          upperBounds,     ... % Parameter upper bounds
    'options',     opts);                % Options defined above

% Create a GlobalSearch object
gs = GlobalSearch;
   
% Run it, returning the best-fitting parameter values and the negative-
% log-likelihood returned by the objective function


[fits(ii,:), nllk] = run(gs,problem);




%write a function to compute the log likelihood

function logL = logLikelihood(adjustedRT, mu, sigma)
    % adjustedRT: Adjusted reaction time (data point)
    % mu: Mean parameter of the LATER model
    % sigma: Standard deviation parameter of the LATER model

    % Compute the probability density function (PDF) of a Gaussian distribution
    pdf = 1 / (sigma * sqrt(2 * pi)) * exp(-(adjustedRT - mu).^2 / (2 * sigma.^2));

    % Compute the log-likelihood
    logL = log(pdf);
end


%final thoughts, cannot be evaluated since the ii function is not defined.
%What is this supposed to represent? According to chatGPT it should be a
%variable defined elsewhere, but I don't see that in any of the code that
%was provided. And I also don't see information regarding where we
%should've used that specific variable in the code we produced ourselves...





