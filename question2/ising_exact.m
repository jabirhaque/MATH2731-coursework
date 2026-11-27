function [probabilities, energies, magnetisation] = ising_exact(N, T)

num_configs = 2^N;

binary_configs = dec2bin(0:num_configs - 1, N) - '0';

configurations = 2 * binary_configs - 1;

energies = zeros(1, num_configs);
magnetisation = zeros(1, num_configs);

for k = 1:num_configs
    spins = configurations(k, :);

    energies(k) = ising_energy(spins);

    magnetisation(k) = ising_magnetisation(spins);
end

weights = exp(-(energies - min(energies))  / T);

Z = sum(weights);


probabilities = weights ./ Z;

end