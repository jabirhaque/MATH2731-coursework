function [magnetisations, energies, final_spins, acceptance_rate] = metropolis_ising(N, T, M, initial_spins)

spins = initial_spins;

magnetisations = zeros(1, M);
energies = zeros(1, M);

accepted = 0;
current_energy = ising_energy(spins);

for k = 1:M
    i = randi(N);
    proposed_spins = spins;
    proposed_spins(i) = -proposed_spins(i);
    
    proposed_energy = ising_energy(proposed_spins);
    delta_energy = proposed_energy - current_energy;
    if delta_energy <=0 || rand < exp(-delta_energy/T)
        spins = proposed_spins;
        accepted = accepted + 1;
        current_energy = proposed_energy;
    end
    magnetisations(k) = ising_magnetisation(spins);
    energies(k) = current_energy;
end

final_spins = spins;
acceptance_rate = accepted/M;

end