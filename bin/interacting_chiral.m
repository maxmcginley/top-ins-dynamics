clc;
clear;

addpath(fullfile(pwd,'..','interacting'));
addpath(fullfile(pwd,'..','TI'));

%********************INPUT DATA****************
sites = 8;
fill = sites/2;
spintot = fill;
U = 1.5;
t_a = 0.9;
t_b = 1.1;
spinful = false;
%********************************************

hoppings = repmat([t_a,t_b],1,sites/2);

f = InteractingTopologicalInsulator(sites,fill,spintot,spinful);

ham = f.interaction_hamiltonian(U) + f.kinetic_hamiltonian(-hoppings);

[Evec,Eval] = eigs(ham,f.dim,'sa');

diag(Eval)

corrmat = eye(sites) - f.correlation_matrix(Evec(:,1));
[k_corrs,ks] = TopologicalInsulator.k_space_matrix(corrmat,2);

k_eigs = zeros(sites,1);

for j = 1:(sites/2)
    k_eigs((2*j-1):(2*j),1) = eig(k_corrs{j});
end

[sort(k_eigs), sort(eig(corrmat))]

invar = TopologicalInsulator.invariant_from_corrmats(k_corrs,ks)