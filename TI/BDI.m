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
sites = 120;
open = false;
hopping_A = 0.3*exp(0i);
hopping_B = 1*exp(0i);
hopping_3A = 0.4*exp(0i);
hopping_3B = 2*exp(0i);
hopping_3A_post = 0.3*exp(-0.5i);
hopping_3B_post = 1.4*exp(0.4i);
times = 0:0.1:5;
k_times = [0:0.01:0.09, 0.1:0.1:10];
site1 = 1;
site2 = 60;
current_site = 20;
cell_size = 2;
hopping_range = 1;
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

cells = sites/cell_size;

ins = TopologicalInsulator_BDI(hopping_A,hopping_B,hopping_3A,hopping_3B,sites,open);
ins2 = TopologicalInsulator_BDI(hopping_A,hopping_B,hopping_3A_post,hopping_3B_post,sites,open);

ins.spectrum
ins2.spectrum

spectra1 = zeros(site2 - site1 + 1,numel(times));
%init_mat = diag(repmat([1,0],1,sites/2));
init_mat = ins.half_filled_correlation_matrix(0);

for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t_1 = ins2.time_evolve_correlation_matrix(init_mat,t);
    spectra1(:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(...
        corrmat_t_1, site1, site2)));
end

%% Test symmetries

[init_trs,init_phs,init_chi] = TopologicalInsulator_SSH.test_symmetries(init_mat);
[t_trs,t_phs,t_chi] = TopologicalInsulator_SSH.test_symmetries(corrmat_t_1);

crit = 1.e-6;

fprintf('Initially: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
fprintf('Time t: PHS = %d ; TRS = %d ; CHI = %d \n',abs(t_phs) < crit,abs(t_trs) < crit,abs(t_chi) < crit);

%% Plotting

figure_handles{end+1} = figure('Name','Entanglement spectrum');
plot(times,spectra1);
ylim([-2,2]);

%figure_handles{end+1} = TopologicalInsulator_SSH.plot_current_entanglement(k_times,top_invars,deriv_top_invars,times,charges,currents,spectra1);

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