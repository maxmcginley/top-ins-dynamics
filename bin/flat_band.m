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
open =  false;
en_scale = 2*pi;
times = 0:0.01:10;
site1 = 5;
site2 = 8;
current_site = 25;
cell_size = 4;
hopping_range = 3;
%*********************************************

if mod(site1,cell_size) ~= 1 || mod(site2,cell_size) ~= 0
    error('Entanglement cut must not be within a unit cell');
end
if mod(sites,cell_size) ~=0 && ~open
    error('Number of sites must be zero modulo cell size for closed system');
end
cells = sites/cell_size;

ham = TopologicalInsulator_Flat(cells,open);

ins_1 = TopologicalInsulator_Flat(en_scale,+1,sites,open);
ins_2 = TopologicalInsulator_Flat(en_scale,-1,sites,open);


zero_wavefunction = ins.zero_mode_wavefunction(-0.01,0.01);
wf_t = ins.time_evolve_creation_operator(zero_wavefunction.',times(end));
zero_wavefunction = [1,zeros(1,sites-1)];

params_out = FlatBandHamiltonian.commensurate_parameters((1/7)*(2*pi),(1/3)*(2*pi),a2,b2,c2,s2);

ham2 = FlatBandHamiltonian.four_cell_hamiltonian(params_out,breaking2,chiral_breaking2,cells,open);
ins2 = TopologicalInsulator(ham2);

evals = eig(ham2)
[t1,t2] = FlatBandHamiltonian.energies_from_parameters(params_out)

t = TopologicalInsulator.validate_current_operator(ham2,site1,site2,hopping_range,open)

init_mat = ins.half_filled_correlation_matrix(-0.0005);
figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(init_mat,open);

occupations = zeros(1,numel(times));
spectra = zeros(site2 - site1 + 1,numel(times));
entropies = zeros(1,numel(times));

currents = zeros(1,numel(times));
% bloch_currents = zeros(1,numel(times));
qubits = zeros(1,numel(times));
wannier_centres = zeros(sites/2,numel(times));

H2 = @(x) -x*exp(-x);

curr_op = TopologicalInsulator.current_operator_from_hamiltonian_test(ham2,current_site,hopping_range,open);

kvals = (2*pi)*(0:(cells-1))/cells;
eig_curr_ks = TopologicalInsulator.k_current_operator(curr_op,kvals,cell_size);

chern_densities = TopologicalInsulator.local_chern_densities(init_mat,cell_size);

prev_wannier = NaN;

for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t = ins2.time_evolve_correlation_matrix(init_mat,t);
    occupations(1,t_index) = TopologicalInsulator.two_particle_expectation(...
        corrmat_t, zero_wavefunction, conj(zero_wavefunction));
    spectra(:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(...
        corrmat_t, site1, site2)));
    entropies(1,t_index) = TopologicalInsulator.entanglement_entropy_from_correlation_matrix(...
        corrmat_t, site1, site2);
    
    %currents(1,t_index) = sum(sum(curr_op .* corrmat_t));
    currents(1,t_index) = trace(curr_op * corrmat_t);
%     bloch_currents(1,t_index) = TopologicalInsulator.bloch_current_vals(ham2,corrmat_t,cell_size,2);
    qubits(1,t_index) = ins2.noneq_qubit_function(init_mat,1,t);
    wannier_centres(:,t_index) = TopologicalInsulator.wannier_centres(corrmat_t,prev_wannier,cell_size);
    prev_wannier = wannier_centres(:,t_index);
end

integrated_current_induced = cumtrapz(times,currents);
% integrated_bloch_current_induced = cumtrapz(times,bloch_currents);
integrated_current_induced(end)
% integrated_bloch_current_induced(end)

[berry_curvatures,kvals] = TopologicalInsulator.berry_curvatures_test(ham,ham2,cell_size,2,1,times);


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
hold on;
plot(times,currents);
% plot(times,bloch_currents);
hold off;
subplot(2,1,2)
hold on;
plot(times,integrated_current_induced);
% plot(times,integrated_bloch_current_induced);
hold off;

figure_handles{end+1} = figure('Name','Qubit preservation');

plot(times,abs(qubits)/abs(qubits(1)));

figure_handles{end+1} = figure('Name','Real space chern density');

plot(1:(sites/cell_size),chern_densities);

figure_handles{end+1} = figure('Name','Correlation matrix zero time');
pos_values = exp(((1:(sites)) + 0.5)*2i*pi/sites);
imagesc(abs(init_mat*diag(pos_values)*init_mat))

figure_handles{end+1} = figure('Name','Wannier centre flow');
plot(times,wannier_centres);
ylim([-5,5]);

figure_handles{end+1} = figure('Name','Current operator in k space');

plot(kvals,eig_curr_ks);

figure_handles{end+1} = figure('Name','Berry Curvatures');
surf(kvals/(2*pi),times,sum(berry_curvatures,3));
view(2);
shading interp;