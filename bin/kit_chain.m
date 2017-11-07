if exist('figure_handles','var') 
    for j = 1:numel(figure_handles)
        if ishandle(figure_handles{j})
            close(figure_handles{j});
        end
    end
    clear('figure_handles');
end

clc;
clear;

figure_handles = cell(1,1);

addpath(fullfile(pwd,'..','TI'));

%******************INPUT DATA*******************
sites = 100;
open = false;
mu_A = 0.2;
t_A = 1*exp(0i);
delta_A = 2i;
mu_B = 2.5;
t_B = 0.5i;
delta_B = 1;
times = 0:0.05:4;
site1 = 27;
site2 = 128;
current_site = 50;
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

ham = TopologicalInsulator.kitaev_chain_hamiltonian(mu_A,t_A,delta_A,sites,open);



ins = TopologicalInsulator(ham);
figure_handles{end+1} = ins.plot_orbitals_in_energy_range(-0.01,0.01);

zero_wavefunction = ins.zero_mode_wavefunction(-0.01,0.01);
wf_t = ins.time_evolve_creation_operator(zero_wavefunction.',times(end));
%zero_wavefunction = [1,zeros(1,2*sites-1)];

ham2 = TopologicalInsulator.kitaev_chain_hamiltonian(mu_B,t_B,delta_B,sites,open);
ins2 = TopologicalInsulator(ham2);

ham_conj = TopologicalInsulator.nambu_operator(1,2*sites) * conj(ham2) * TopologicalInsulator.nambu_operator(1,2*sites);
sum(sum(abs(ham2 + ham_conj)))

init_mat = ins.half_filled_correlation_matrix(-0.01);
figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(init_mat,open);

occupations = zeros(1,numel(times));
spectra = zeros(site2 - site1 + 1,numel(times));
entropies = zeros(1,numel(times));
currents = zeros(1,numel(times));

H2 = @(x) -x*exp(-x);

curr_op = TopologicalInsulator.current_operator_from_hamiltonian(ham2,current_site,open);

for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t = ins2.time_evolve_correlation_matrix(init_mat,t);
    occupations(1,t_index) = TopologicalInsulator.two_particle_expectation(...
        corrmat_t, zero_wavefunction, conj(zero_wavefunction));
    spectra(:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(...
        corrmat_t, site1, site2)));
    entropies(1,t_index) = TopologicalInsulator.entanglement_entropy_from_correlation_matrix(...
        corrmat_t, site1, site2);
    
    currents(1,t_index) = sum(sum(curr_op .* corrmat_t));
end

integrated_current_induced = cumtrapz(currents,ones(1,numel(times)).*(times(2) - times(1)));

figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(corrmat_t,open);

topological_spectrum_indices = [(site2 - site1 + 1)/2, (site2 - site1 + 3)/2];

%% Testing symmetry

[init_trs,init_phs,init_chi] = TopologicalInsulator.test_symmetries(init_mat);
[t_trs,t_phs,t_chi] = TopologicalInsulator.test_symmetries(corrmat_t);

crit = 1.e-6;

fprintf('Initially: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
fprintf('Time t: PHS = %d ; TRS = %d ; CHI = %d \n',abs(t_phs) < crit,abs(t_trs) < crit,abs(t_chi) < crit);

%% Time evolved ham

figure_handles{end+1} = figure('Name','Hamiltonian zero time');
imagesc(abs(ham))
figure_handles{end+1} = figure('Name','Hamiltonian inter time');
ham_int = ins2.time_evolve_hamiltonian(ham,times(uint32(floor(numel(times)/2))));
imagesc(abs(ham_int))
figure_handles{end+1} = figure('Name','Hamiltonian full time');
ham_t = ins2.time_evolve_hamiltonian(ham,times(end));
imagesc(abs(ham_t))

%% Plotting

figure_handles{end+1} = figure('Name','Occupation evolution');

plot(times,real(occupations));

figure_handles{end+1} = figure('Name','Entanglement spectra evolution');

plot(times,sort(real(spectra(topological_spectrum_indices,:))));

figure_handles{end+1} = figure('Name','Entanglement entropy');

plot(times,entropies);

figure_handles{end+1} = figure('Name','Time evolved zero wavefunction');

plot(1:(2*sites),abs(wf_t));

figure_handles{end+1} = figure('Name','Current induced');

subplot(2,1,1);
plot(times,currents);
subplot(2,1,2)
plot(times,integrated_current_induced);