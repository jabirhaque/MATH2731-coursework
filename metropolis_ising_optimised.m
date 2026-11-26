function [magnetisations, energies, final_spins, acceptance_rate] = metropolis_ising(N, T, M, initial_spins)

spins = initial_spins;

magnetisations = zeros(1, M);
energies = zeros(1, M);

accepted = 0;
current_energy = ising_energy(spins);

current_spin_sum = sum(spins);

for k = 1:M
    i = randi(N);
    delta_energy = ising_delta_energy(spins, i);
    if delta_energy <=0 || rand < exp(-delta_energy/T)
        current_spin_sum = current_spin_sum - 2 * spins(i);
        spins(i) = -spins(i);
        accepted = accepted + 1;
        current_energy = current_energy + delta_energy;
    end
    magnetisations(k) = current_spin_sum / N;
    energies(k) = current_energy;
end

final_spins = spins;
acceptance_rate = accepted/M;

end