if exist('figure_handles','var') 
    for j = 1:numel(figure_handles)
        if ishandle(figure_handles{j})
            close(figure_handles{j});
        end
    end
    clear('figure_handles');
end

%%
% $x^2+e^{\pi i}$
%%
% 
% $$e^{\pi i} + 1 = 0$$
% 

clc;
clear;

figure_handles = cell(1,1);

addpath(fullfile(pwd,'..','TI'));

%******************INPUT DATA*******************
sites = 24;
open = true;
hopping_A = 0.3*exp(0.0i);
hopping_B = 1*exp(0.4i);
hopping_A2 = 0.8*exp(0.4i);
hopping_B2 = 1*exp(0.0i);
hopping_A3 = 0.3;
hopping_B3 = 1;
hopping_A4 = 0.8;
hopping_B4 = 1;
hopping_A_third = 0.05;
hopping_B_third = 0.05;
times = [0:0.01:0.09, 0.1:0.1:10];
site1 = 7;
site2 = 18;
current_site = 12;
cell_size = 2;
hopping_range = 3;
haldane_data_directory = 'haldane_data';
%*********************************************

header_lines = 22;

haldane_data_1 = importdata(fullfile(haldane_data_directory,'data_trs_axial.txt'),',',header_lines);
haldane_data_2 = importdata(fullfile(haldane_data_directory,'data_trs_only.txt'),',',header_lines);

h_times_1 = haldane_data_1.data(:,1).';
h_times_2 = haldane_data_2.data(:,1).';
h_ents_1 = -log(haldane_data_1.data(:,3:end)).';
h_ents_2 = -log(haldane_data_2.data(:,3:end)).';

hal_times = 0:0.01:2;

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

cells = sites/cell_size;

ins1 = TopologicalInsulator_SSH(hopping_A,hopping_B,sites,open,hopping_A_third,hopping_B_third);
ins1.spectrum
ins2 = TopologicalInsulator_SSH(hopping_A2,hopping_B2,sites,open,hopping_A_third,hopping_B_third);
ins3 = TopologicalInsulator_SSH(hopping_A3,hopping_B3,sites,open,hopping_A_third,hopping_B_third);
ins4 = TopologicalInsulator_SSH(hopping_A4,hopping_B4,sites,open,hopping_A_third,hopping_B_third);

occupations = zeros(1,numel(times));

currents_2 = zeros(1,numel(times));
charges_2 = zeros(1,numel(times));
currents_4 = zeros(1,numel(times));
charges_4 = zeros(1,numel(times));

spec_1 = zeros(1,numel(times));
spec_2 = zeros(1,numel(times));

curr_op_2 = TopologicalInsulator.current_operator_from_hamiltonian(ins2.hamiltonian,current_site,hopping_range);
curr_op_4 = TopologicalInsulator.current_operator_from_hamiltonian(ins4.hamiltonian,current_site,hopping_range);


% charge_op = diag([zeros(1,current_site),ones(1,sites-current_site)]);
charge_op = diag([ones(1,current_site),zeros(1,sites-current_site)]);

%init_mat = diag(repmat([1,0],1,sites/2));
init_mat_1 = ins1.half_filled_correlation_matrix(-0.1);
init_mat_3 = ins3.half_filled_correlation_matrix(-0.1);

init_charge_1 = real(trace(charge_op * init_mat_1.')); init_charge_1 = floor(init_charge_1);
init_charge_3 = real(trace(charge_op * init_mat_3.')); init_charge_3 = floor(init_charge_3);


for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t_2 = ins2.time_evolve_correlation_matrix(init_mat_1,t);
    corrmat_t_4 = ins4.time_evolve_correlation_matrix(init_mat_3,t);
    currents_2(1,t_index) = real(trace(curr_op_2 * corrmat_t_2.'));
    charges_2(1,t_index) = real(trace(charge_op * corrmat_t_2.')) - init_charge_1;
    
    currents_4(1,t_index) = real(trace(curr_op_4 * corrmat_t_4.'));
    charges_4(1,t_index) = real(trace(charge_op * corrmat_t_4.')) - init_charge_1;
    
    spec_1(1,t_index) = min(abs(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(corrmat_t_2, site1, site2)));
    spec_2(1,t_index) = min(abs(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(corrmat_t_4, site1, site2)));
end

integrated_current_induced = (times(2) - times(1))*cumtrapz(currents_2);

topological_spectrum_indices = [(site2 - site1 + 1)/2, (site2 - site1 + 3)/2];

num_ks = 100;
k_vals = 2*pi * (0:(num_ks - 1))/ num_ks;
init_spinors_1 = ins1.BL_ground_spinors(k_vals);
init_spinors_3 = ins3.BL_ground_spinors(k_vals);
[top_invars_1] = ins2.BL_topological_invariant(init_spinors_1,times,k_vals);
[top_invars_2] = ins4.BL_topological_invariant(init_spinors_3,times,k_vals);


%% Test symmetries

[init_trs,init_phs,init_chi] = TopologicalInsulator_SSH.test_symmetries(eye(size(init_mat_1)) - 2*init_mat_1);
[t_trs,t_phs,t_chi] = TopologicalInsulator_SSH.test_symmetries(corrmat_t_2);

crit = 1.e-6;

fprintf('Initially: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
fprintf('Time t: PHS = %d ; TRS = %d ; CHI = %d \n',abs(t_phs) < crit,abs(t_trs) < crit,abs(t_chi) < crit);

%% Plotting

figure_handles{end+1} = TopologicalInsulator_SSH.plot_current_entanglement(times,top_invars_1,top_invars_2,charges_2,charges_4,spec_1,spec_2,h_times_1,h_ents_1,h_times_2,h_ents_2);

% 
% figure_handles{end+1} = TopologicalInsulator.plot_entanglement_spectrum(spectra1,spectra2,times);
% 
% figure_handles{end+1} = figure('Name','Current induced');
% subplot(2,1,1);
% hold on;
% plot(times,currents);
% plot(times,currents_a);
% hold off;
% subplot(2,1,2)
% plot(times,integrated_current_induced);
% 
% figure_handles{end+1} = figure('Name','Topological Invariant');
% hold on;
% yyaxis left;
% plot(times,top_invars);
% yyaxis right;
% plot(times(1:(end-1)),-deriv_top_invars);
% hold off;
% xlabel('Time $t$','interpreter','latex');
% ylabel('$\nu$','interpreter','latex');

%% Animation

% figure_handles{end+1} = chern_density_animation(0:(cells-1),chern_densities,times,5);

% figure_handles{end+1} = figure('Name','Wannier centre flow');
% plot(times,wannier_centres);
% ylim([-5,5]);

% figure_handles{end+1} = figure('Name','Berry Curvatures');
% surf(kvals/(2*pi),times(1:end-1),berry_curvatures);
% view(2);
% shading interp;