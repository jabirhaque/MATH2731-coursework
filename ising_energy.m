function E = ising_energy(spins)
E = -sum(spins .* circshift(spins, -1));
end