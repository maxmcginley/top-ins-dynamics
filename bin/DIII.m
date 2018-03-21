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
sites = 40;
open = true;
mu = 0.5;
delp = 1;
dels_1 = 0.2;
dels_2 = -0.8;
alpha = 0.3;
times = 0:0.1:5;
k_times = [0:0.01:0.09, 0.1:0.1:10];
site1 = 1;
site2 = 20;
current_site = 20;
cell_size = 2;
hopping_range = 1;
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

cells = sites/cell_size;

ins = TopologicalInsulator_DIII(mu,delp,dels_1,alpha,sites,open);
ins2 = TopologicalInsulator_DIII(mu,delp,dels_2,alpha,sites,open);

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