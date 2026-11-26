00000000000%% Metropolis-Hastings Sampling on [0, 1]
% Target: p(x) proportional to exp(-5*(x - 0.5)^2) * cos(3*pi*x) + 1 on [0, 1]

clear; clc; close all;
rng(42); % For reproducibility

%% 1. Define Unnormalized Target Density and Normalizing Constant
p_tilde = @(x) (x >= 0 & x <= 1) .* (exp(-5 * (x - 0.5).^2) .* cos(3 * pi * x) + 1);

% Compute normalizing constant Z via numerical integration (for plotting true density)
Z = integral(p_tilde, 0, 1);
p_true = @(x) p_tilde(x) / Z;

%% 2. Run Metropolis-Hastings Sampler
N = 50000;          % Total number of iterations
burn_in = 5000;     % Burn-in samples to discard
sigma = 0.2;        % Proposal standard deviation (Gaussian random walk)
x0 = 0.05;          % Initial state near boundary to observe burn-in/exploration

samples = zeros(N, 1);
samples(1) = x0;
accepted = 0;

for t = 1:(N - 1)
    x_curr = samples(t);
    
    % Symmetric Gaussian proposal: x* ~ N(x_curr, sigma^2)
    x_prop = x_curr + sigma * randn();
    
    % Acceptance ratio (symmetric proposal cancels q terms; out-of-bounds has p_tilde = 0)
    if x_prop >= 0 && x_prop <= 1
        alpha = min(1, p_tilde(x_prop) / p_tilde(x_curr));
    else
        alpha = 0;
    end
    
    % Accept or reject
    if rand() <= alpha
        samples(t + 1) = x_prop;
        accepted = accepted + 1;
    else
        samples(t + 1) = x_curr;
    end
end

accept_rate = accepted / (N - 1);
fprintf('Acceptance rate (sigma = %.2f): %.2f%%\n', sigma, 100 * accept_rate);

%% 3. Visualize Chain Exploration & Compare Histogram to True Density
figure('Position', [100, 100, 1000, 750]);

% Subplot 1: Trace plot of the first 1,000 steps (shows exploration across [0, 1])
subplot(2, 2, 1);
plot(1:1000, samples(1:1000), 'LineWidth', 0.8);
ylim([0, 1]);
xlabel('Iteration');
ylabel('x_t');
title('Trace Plot (First 1,000 Iterations)');
grid on;

% Subplot 2: Full trace plot
subplot(2, 2, 2);
plot(1:N, samples, 'LineWidth', 0.3);
ylim([0, 1]);
xlabel('Iteration');
ylabel('x_t');
title(sprintf('Full Trace Plot (Acceptance Rate = %.1f%%)', 100 * accept_rate));
grid on;

% Subplot 3: Normalized histogram vs. true density (post burn-in)
subplot(2, 2, [3, 4]);
post_burnin = samples((burn_in + 1):end);
histogram(post_burnin, 50, 'Normalization', 'pdf', ...
    'FaceColor', [0.3, 0.6, 0.85], 'EdgeColor', 'white');
hold on;

x_grid = linspace(0, 1, 500);
plot(x_grid, p_true(x_grid), 'r-', 'LineWidth', 2.2);
hold off;

xlim([0, 1]);
xlabel('x');
ylabel('Probability Density');
title(sprintf('Sample Histogram (N = %d, Burn-in = %d) vs. True Density p(x)', N, burn_in));
legend('MH Sample Histogram', 'True Normalized Density p(x)', 'Location', 'northeast');
grid on;