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
sites = 50;
open = false;
hopping_A = 0.4;
hopping_B = 1*exp(0i);
hopping_A2 = 0.4*exp(0.5i);
hopping_B2 = 2.5*exp(-0.2i);
times = 0:0.5:40;
site1 = 3;
site2 = 8;
current_site = 5;
cell_size = 2;
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

ham = TopologicalInsulator.SSH_hamiltonian(hopping_A,hopping_B,sites,open);

ins = TopologicalInsulator(ham);
figure_handles{end+1} = ins.plot_orbitals_in_energy_range(-0.01,0.01);

zero_wavefunction = ins.zero_mode_wavefunction(-0.01,0.01);
wf_t = ins.time_evolve_creation_operator(zero_wavefunction.',times(end));
zero_wavefunction = [1,zeros(1,sites-1)];

ham2 = TopologicalInsulator.SSH_hamiltonian(hopping_A2,hopping_B2,sites,open);
ins2 = TopologicalInsulator(ham2);



init_mat = ins.half_filled_correlation_matrix(0.01);
figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(init_mat,open);

chern_densities = TopologicalInsulator.local_chern_densities(init_mat,cell_size);

occupations = zeros(1,numel(times));
spectra = zeros(site2 - site1 + 1,numel(times));
entropies = zeros(1,numel(times));

currents = zeros(1,numel(times));
qubits = zeros(1,numel(times));
wannier_centres = zeros(sites/2,numel(times));
prev_wannier = NaN;

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
    wannier_centres(:,t_index) = TopologicalInsulator.wannier_centres(corrmat_t,prev_wannier,cell_size);

    qubits(1,t_index) = ins2.noneq_qubit_function(init_mat,1,t);
end

integrated_current_induced = (times(2) - times(1))*cumtrapz(currents);

figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(corrmat_t,open);

topological_spectrum_indices = [(site2 - site1 + 1)/2, (site2 - site1 + 3)/2];

%% Test symmetries

[init_trs,init_phs,init_chi] = TopologicalInsulator.test_symmetries(init_mat);
[t_trs,t_phs,t_chi] = TopologicalInsulator.test_symmetries(corrmat_t);

crit = 1.e-6;

fprintf('Initially: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
fprintf('Time t: PHS = %d ; TRS = %d ; CHI = %d \n',abs(t_phs) < crit,abs(t_trs) < crit,abs(t_chi) < crit);

%% Plotting

figure_handles{end+1} = figure('Name','Occupation evolution');

plot(times,real(occupations));

figure_handles{end+1} = figure('Name','Entanglement spectra evolution');

plot(times,sort(real(spectra(topological_spectrum_indices,:))));

figure_handles{end+1} = figure('Name','Entanglement entropy');

plot(times,entropies);

figure_handles{end+1} = figure('Name','Time evolved zero wavefunction');

plot(1:(sites),abs(wf_t));

figure_handles{end+1} = figure('Name','Current induced');

subplot(2,1,1);
plot(times,currents);
subplot(2,1,2)
plot(times,integrated_current_induced);

figure_handles{end+1} = figure('Name','Qubit preservation');

plot(times,abs(qubits)/abs(qubits(1)));

figure_handles{end+1} = figure('Name','Real space chern density');

plot(1:(sites/cell_size),chern_densities);

figure_handles{end+1} = figure('Name','Wannier centre flow');
plot(times,wannier_centres);
ylim([-5,5]);