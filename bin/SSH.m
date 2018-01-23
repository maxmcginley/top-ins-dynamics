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
hopping_A = 0.3*exp(0i);
hopping_B = 1*exp(0i);
hopping_A2 = 0.7*exp(-0.1i);
hopping_B2 = 1*exp(0.1i);
hopping_A3 = 0.7*exp(-0i);
hopping_B3 = 1*exp(0i);
times = 0:0.1:20;
site1 = 25;
site2 = 74;
current_site = 49;
cell_size = 2;
hopping_range = 1;
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

cells = sites/cell_size;

ham = TopologicalInsulator.SSH_hamiltonian(hopping_A,hopping_B,sites,open);
ins = TopologicalInsulator(ham);
figure_handles{end+1} = ins.plot_orbitals_in_energy_range(-0.01,0.01);

zero_wavefunction = ins.zero_mode_wavefunction(-0.01,0.01);
wf_t = ins.time_evolve_creation_operator(zero_wavefunction.',times(end));
zero_wavefunction = [1,zeros(1,sites-1)];

ham2 = TopologicalInsulator.SSH_hamiltonian(hopping_A2,hopping_B2,sites,open);
ins2 = TopologicalInsulator(ham2);

ham3 = TopologicalInsulator.SSH_hamiltonian(hopping_A3,hopping_B3,sites,open);
ins3 = TopologicalInsulator(ham3);



init_mat = ins.half_filled_correlation_matrix(0.01);
figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(init_mat,open);


occupations = zeros(1,numel(times));
spectra1 = zeros(site2 - site1 + 1,numel(times));
spectra2 = zeros(site2 - site1 + 1,numel(times));
entropies = zeros(1,numel(times));
chern_densities = zeros(numel(times),cells);

currents = zeros(1,numel(times));
qubits = zeros(1,numel(times));
wannier_centres = zeros(sites/2,numel(times));
prev_wannier = NaN;

H2 = @(x) -x*exp(-x);

curr_op = TopologicalInsulator.current_operator_from_hamiltonian_test(ham2,current_site,hopping_range,open);

for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t_1 = ins2.time_evolve_correlation_matrix(init_mat,t);
    corrmat_t_2 = ins3.time_evolve_correlation_matrix(init_mat,t);
%     occupations(1,t_index) = TopologicalInsulator.two_particle_expectation(...
%         corrmat_t, zero_wavefunction, conj(zero_wavefunction));
    spectra1(:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(...
        corrmat_t_1, site1, site2)));
    spectra2(:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(...
        corrmat_t_2, site1, site2)));
%     entropies(1,t_index) = TopologicalInsulator.entanglement_entropy_from_correlation_matrix(...
%         corrmat_t, site1, site2);
%     chern_densities(t_index,:) = TopologicalInsulator.local_chern_densities(corrmat_t,cell_size);
    
    currents(1,t_index) = trace(curr_op * corrmat_t_1);
    wannier_centres(:,t_index) = TopologicalInsulator.wannier_centres(corrmat_t_1,prev_wannier,cell_size);

%     qubits(1,t_index) = ins2.noneq_qubit_function(init_mat,1,t);
end

% ham_b = TopologicalInsulator.SSH_hamiltonian(hopping_A,hopping_B,1000,false);
% ham2_b = TopologicalInsulator.SSH_hamiltonian(hopping_A2,hopping_B2,1000,false);
%[berry_curvatures,kvals] = TopologicalInsulator.berry_curvatures(ham_b,ham2_b,cell_size,1,1,times);


integrated_current_induced = (times(2) - times(1))*cumtrapz(currents);

figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(corrmat_t_1,open);

topological_spectrum_indices = [(site2 - site1 + 1)/2, (site2 - site1 + 3)/2];

%% Test symmetries

[init_trs,init_phs,init_chi] = TopologicalInsulator.test_symmetries(init_mat);
[t_trs,t_phs,t_chi] = TopologicalInsulator.test_symmetries(corrmat_t_1);

crit = 1.e-6;

fprintf('Initially: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
fprintf('Time t: PHS = %d ; TRS = %d ; CHI = %d \n',abs(t_phs) < crit,abs(t_trs) < crit,abs(t_chi) < crit);

%% Plotting

figure_handles{end+1} = figure('Name','Occupation evolution');

plot(times,real(occupations));

figure_handles{end+1} = TopologicalInsulator.plot_entanglement_spectrum(spectra1,spectra2,times);

% figure_handles{end+1} = figure('Name','Time evolved zero wavefunction');
% 
% plot(1:(sites),abs(wf_t));
% 
figure_handles{end+1} = figure('Name','Current induced');

subplot(2,1,1);
plot(times,currents);
subplot(2,1,2)
plot(times,integrated_current_induced);
% 
% figure_handles{end+1} = figure('Name','Qubit preservation');
% 
% plot(times,abs(qubits)/abs(qubits(1)));

%% Animation

% figure_handles{end+1} = chern_density_animation(0:(cells-1),chern_densities,times,5);

figure_handles{end+1} = figure('Name','Wannier centre flow');
plot(times,wannier_centres);
ylim([-5,5]);

% figure_handles{end+1} = figure('Name','Berry Curvatures');
% surf(kvals/(2*pi),times(1:end-1),berry_curvatures);
% view(2);
% shading interp;