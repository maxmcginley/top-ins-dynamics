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
h01 = [0;0;1];
hc1 = [0;0;0];
hs1 = [0;0;0];
%[hc1,hs1]= FlatBandHamiltonian.cos_sin_vectors(h01,hc1);
breaking1 = 5;
chiral_breaking1 = 0;

h02 = [0;0;1];
hc2 = [1;0;0];
hs2 = [0;-1;0];
%[hc2,hs2]= FlatBandHamiltonian.cos_sin_vectors(h02,hc2);
period = 6;
breaking2 = 0;
chiral_breaking2 = 0;
times = 0:0.02:12;
site1 = 5;
site2 = 8;
current_site = 11;
cell_size = 2;
hopping_range = 3;
%*********************************************

cells = sites/cell_size;

if mod(site1,cell_size) ~= 1 || mod(site2,cell_size) ~= 0
    error('Entanglement cut must not be within a unit cell');
end
if mod(sites,cell_size) ~=0 && ~open
    error('Number of sites must be zero modulo cell size for closed system');
end
cells = sites/cell_size;

ham = FlatBandHamiltonian.two_cell_hamiltonian(h01,hc1,hs1,cells,open);
eig(ham)

ins = TopologicalInsulator(ham);

ham2 = FlatBandHamiltonian.two_cell_hamiltonian(h02,hc2,hs2,cells,open);


evals = eig(ham2);
ham2 = ham2 * 2*pi/(period*abs(evals(1)));

ins2 = TopologicalInsulator(ham2);

init_mat = ins.half_filled_correlation_matrix(0.01);
%figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(init_mat,open);

currents = zeros(1,numel(times));
% bloch_currents = zeros(1,numel(times));
wannier_centres = zeros(sites/2,numel(times));

H2 = @(x) -x*exp(-x);

curr_op = TopologicalInsulator.current_operator_from_hamiltonian_test(ham2,current_site,hopping_range,open);

kvals = (2*pi)*(0:(cells-1))/cells;
eig_curr_ks = TopologicalInsulator.k_current_operator(curr_op,kvals,cell_size);

%chern_densities = TopologicalInsulator.local_chern_densities(init_mat,cell_size);

prev_wannier = NaN;

for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t = ins2.time_evolve_correlation_matrix(init_mat,t);
    
    currents(1,t_index) = trace(curr_op * corrmat_t);
    %wannier_centres(:,t_index) = TopologicalInsulator.wannier_centres(corrmat_t,prev_wannier,cell_size);
    %prev_wannier = wannier_centres(:,t_index);
end

[berry_curvatures,kvals] = TopologicalInsulator.berry_curvatures_test(ham,ham2,cell_size,1,1,times);

integrated_current_induced = cumtrapz(times,currents);
% integrated_bloch_current_induced = cumtrapz(times,bloch_currents);
integrated_current_induced(end)
% integrated_bloch_current_induced(end)


figure_handles{end+1} = TopologicalInsulator.correlation_length_plot(corrmat_t,open);

topological_spectrum_indices = [(site2 - site1 + 1)/2, (site2 - site1 + 3)/2];
% 
% %% Test symmetries
% 
% [init_trs,init_phs,init_chi] = TopologicalInsulator.test_symmetries(init_mat);
% [t_trs,t_phs,t_chi] = TopologicalInsulator.test_symmetries(corrmat_t);
% 
% crit = 1.e-6;
% 
% fprintf('Initially: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
% fprintf('Time t: PHS = %d ; TRS = %d ; CHI = %d \n',abs(t_phs) < crit,abs(t_trs) < crit,abs(t_chi) < crit);

%% Plotting

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


figure_handles{end+1} = figure('Name','Correlation matrix zero time');
pos_values = exp(((1:(sites)) + 0.5)*2i*pi/sites);
imagesc(abs(init_mat*diag(pos_values)*init_mat))

% figure_handles{end+1} = figure('Name','Wannier centre flow');
% plot(times,wannier_centres);
% ylim([-5,5]);

figure_handles{end+1} = figure('Name','Berry Curvatures');
surf(kvals/(2*pi),times,berry_curvatures);
view(2);
shading interp;