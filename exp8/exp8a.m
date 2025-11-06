clc; clear; close all;

% Given Cost Coefficients
a = [500 400 200];      % Constant coefficients
b = [5.3 5.5 5.8];      % Linear coefficients
c = [0.004 0.006 0.009];% Quadratic coefficients

% Total Power Demand (MW)
Pd = 800;

% Initial guess for lambda
lambda = 6.0;

% Convergence tolerance
tolerance = 1e-5;

% Maximum iterations
max_iter = 1000;

% Iteration counter
iter = 0;

while true
    iter = iter + 1;
    
    % Compute power output for each unit from lambda
    P1 = (lambda - b(1)) / (2 * c(1));
    P2 = (lambda - b(2)) / (2 * c(2));
    P3 = (lambda - b(3)) / (2 * c(3));
    
    % Total generated power
    P_total = P1 + P2 + P3;
    
    % Check for convergence (power balance)
    if abs(P_total - Pd) < tolerance
        break;
    end
    
    % Compute dP/dλ (sensitivity)
    dP_dLambda = (1 / (2 * c(1))) + (1 / (2 * c(2))) + (1 / (2 * c(3)));
    
    % Update lambda properly
    lambda = lambda + (Pd - P_total) / dP_dLambda;
    
    if iter > max_iter
        disp('Lambda iteration did not converge.');
        break;
    end
end

% Display results
disp('-----------------------------------------');
disp('      Economic Dispatch Results');
disp('-----------------------------------------');
fprintf('Lambda (Final): %.4f\n', lambda);
fprintf('Generator 1 Output (P1): %.4f MW\n', P1);
fprintf('Generator 2 Output (P2): %.4f MW\n', P2);
fprintf('Generator 3 Output (P3): %.4f MW\n', P3);
fprintf('Total Generation: %.4f MW\n', P_total);
fprintf('Total Demand: %.4f MW\n', Pd);

% Compute total cost
C1 = a(1) + b(1)*P1 + c(1)*P1^2;
C2 = a(2) + b(2)*P2 + c(2)*P2^2;
C3 = a(3) + b(3)*P3 + c(3)*P3^2;
Total_Cost = C1 + C2 + C3;

fprintf('Total Generation Cost: र%.2f\n', Total_Cost);
fprintf('Number of Iterations: %d\n', iter);
disp('-----------------------------------------');
