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
sites = 60;
open = false;
hopping_A = 0.3*exp(0i);
hopping_B = 1*exp(0i);
hopping_A2 = 0.7*exp(-0.1i);
hopping_B2 = 1*exp(0.1i);
hopping_A3 = 0.7*exp(-0i);
hopping_B3 = 1*exp(0i);
times = 0:0.1:20;
site1 = 13;
site2 = 40;
current_site = 18;
cell_size = 2;
hopping_range = 1;
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

cells = sites/cell_size;

ins = TopologicalInsulator_SSH(hopping_A,hopping_B,sites,open);
ins2 = TopologicalInsulator_SSH(hopping_A2,hopping_B2,sites,open);
ins3 = TopologicalInsulator_SSH(hopping_A3,hopping_B3,sites,open);

occupations = zeros(1,numel(times));
spectra1 = zeros(site2 - site1 + 1,numel(times));
spectra2 = zeros(site2 - site1 + 1,numel(times));

currents = zeros(1,numel(times));
currents_a = zeros(1,numel(times));
%wannier_centres = zeros(sites/2,numel(times));
prev_wannier = NaN;

H2 = @(x) -x*exp(-x);

curr_op_2 = TopologicalInsulator.current_operator_from_hamiltonian(ins2.hamiltonian,current_site,hopping_range);
curr_op_2a = TopologicalInsulator.current_operator_from_hamiltonian(ins2.hamiltonian,current_site+4,hopping_range);

%init_mat = ins.half_filled_correlation_matrix(-0.0005);
init_mat = diag(repmat([1,0],1,sites/2));

for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t_1 = ins2.time_evolve_correlation_matrix(init_mat,t);
    corrmat_t_2 = ins3.time_evolve_correlation_matrix(init_mat,t);
    spectra1(:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(...
        corrmat_t_1, site1, site2)));
    spectra2(:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(...
        corrmat_t_2, site1, site2)));
    currents(1,t_index) = trace(curr_op_2 * corrmat_t_1.');
    currents_a(1,t_index) = sum(sum(curr_op_2a .* corrmat_t_1));
    %wannier_centres(:,t_index) = TopologicalInsulator.wannier_centres(corrmat_t_1,prev_wannier,cell_size);
end

integrated_current_induced = (times(2) - times(1))*cumtrapz(currents);

figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(corrmat_t_1,open);

topological_spectrum_indices = [(site2 - site1 + 1)/2, (site2 - site1 + 3)/2];

num_ks = 300;
k_vals = 2*pi * (0:(num_ks - 1))/ num_ks;
[top_invars,k_currs] = ins2.BL_topological_invariant(TopologicalInsulator_SSH.BL_constant_spinor([1;0],k_vals),times,k_vals);
deriv_top_invars = diff(top_invars)./ diff(times);


%% Test symmetries

[init_trs,init_phs,init_chi] = TopologicalInsulator_SSH.test_symmetries(init_mat);
[t_trs,t_phs,t_chi] = TopologicalInsulator_SSH.test_symmetries(corrmat_t_1);

crit = 1.e-6;

fprintf('Initially: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
fprintf('Time t: PHS = %d ; TRS = %d ; CHI = %d \n',abs(t_phs) < crit,abs(t_trs) < crit,abs(t_chi) < crit);

%% Plotting

figure_handles{end+1} = figure('Name','Occupation evolution');
plot(times,real(occupations));

figure_handles{end+1} = TopologicalInsulator.plot_entanglement_spectrum(spectra1,spectra2,times);

figure_handles{end+1} = figure('Name','Current induced');
subplot(2,1,1);
hold on;
plot(times,currents);
plot(times,currents_a);
plot(times,k_currs);
hold off;
subplot(2,1,2)
plot(times,integrated_current_induced);

figure_handles{end+1} = figure('Name','Topological Invariant');
hold on;
yyaxis left;
plot(times,top_invars);
yyaxis right;
plot(times(1:(end-1)),-deriv_top_invars);
plot(times,k_currs);
hold off;
xlabel('Time $t$','interpreter','latex');
ylabel('$\nu$','interpreter','latex');

%% Animation

% figure_handles{end+1} = chern_density_animation(0:(cells-1),chern_densities,times,5);

% figure_handles{end+1} = figure('Name','Wannier centre flow');
% plot(times,wannier_centres);
% ylim([-5,5]);

% figure_handles{end+1} = figure('Name','Berry Curvatures');
% surf(kvals/(2*pi),times(1:end-1),berry_curvatures);
% view(2);
% shading interp;