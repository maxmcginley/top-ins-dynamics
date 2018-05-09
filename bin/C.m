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
sites = 128;
open = false;
mu = 1.2;
t = 1;
del_1 = 0.5;
del_2 = -0.3;
times = 0:0.1:5;
k_times = [0:0.01:0.09, 0.1:0.1:10];
site1 = 1;
site2 = 64;
current_site = 20;
cell_size = 4;
hopping_range = 1;
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

cells = sites/cell_size;

num_ks = 100;
k_vals = 2*pi * (0:(num_ks - 1))/ num_ks;

ins = TopologicalInsulator_C(t,mu,del_1,sites,open);
ins2 = TopologicalInsulator_C(t,mu,del_2,sites,open);


init_spinors = ins.BL_ground_state_spinors(k_vals);
[top_invars_1] = ins2.BL_topological_invariant(init_spinors,times,k_vals);

%% Plotting

figure_handles{end+1} = figure('Name','Entanglement spectrum');
plot(times,top_invars_1);
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