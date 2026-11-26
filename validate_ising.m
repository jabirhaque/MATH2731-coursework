clear;
clc;
rng(42);

N = 10;
T = 2;
M = 100000;
burn_in = 10000;

[probabilities, exact_energies, exact_magnetisations] = ...
    ising_exact(N, T);

exact_mean_energy = sum(probabilities .* exact_energies);
exact_mean_magnetisation = sum(probabilities .* exact_magnetisations);
exact_mean_abs_magnetisation = sum(probabilities .* abs(exact_magnetisations));

initial_spins = 2 * randi([0, 1], 1, N) - 1;

[MH_magnetisations, MH_energies, ~, acceptance_rate] = ...
    metropolis_ising_optimised(N, T, M, initial_spins);

MH_magnetisations = MH_magnetisations(burn_in+1:end);
MH_energies = MH_energies(burn_in+1:end);

MH_mean_energy = mean(MH_energies);

MH_mean_magnetisation = mean(MH_magnetisations);

MH_mean_abs_magnetisation = mean(abs(MH_magnetisations));

energy_error = MH_mean_energy - exact_mean_energy;
magnetisation_error = MH_mean_magnetisation - exact_mean_magnetisation;
abs_magnetisation_error = MH_mean_abs_magnetisation - exact_mean_abs_magnetisation;

quantities = ["Energy"; "Magnetisation"; "Absolute magnetisation"];

exact_values = [
    exact_mean_energy;
    exact_mean_magnetisation;
    exact_mean_abs_magnetisation
    ];

MH_values = [
    MH_mean_energy;
    MH_mean_magnetisation;
    MH_mean_abs_magnetisation
    ];

errors = abs(exact_values - MH_values);

results = table(quantities, exact_values, MH_values, errors, ...
    'VariableNames', {'Quantity', 'Exact', 'MH', 'AbsoluteError'});

disp(results);

mu_values = (2 * (0:N) - N) / N;

exact_mu_probs = zeros(size(mu_values));
MH_mu_probs = zeros(size(mu_values));

for k = 1:length(mu_values)

    exact_mask = abs(exact_magnetisations - mu_values(k)) < 1e-10;

    exact_mu_probs(k) = sum(probabilities(exact_mask));

    MH_mask = abs(MH_magnetisations - mu_values(k)) < 1e-10;

    MH_mu_probs(k) = mean(MH_mask);

end

figure;

bar(mu_values, [exact_mu_probs; MH_mu_probs]');

xlabel('Magnetisation');
ylabel('Probability');
title('Exact vs MH magnetisation distribution');

legend('Exact Boltzmann', 'Metropolis-Hastings');
grid on;

TV_distance = 1/2 * sum(abs(exact_mu_probs - MH_mu_probs));

fprintf('Total variation distance: %.6f\n', TV_distance);