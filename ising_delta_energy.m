function dE = ising_delta_energy(spins, i)
N = length(spins);

left = mod(i-2, N) + 1;
right = mod(i, N) + 1;

dE = 2 * spins(i) * (spins(left) + spins(right));

end