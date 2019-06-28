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
sites = 32;
open = false;
t = 1; %Average hopping matrix element
laser_detuning = 0.2; %Difference between ladder hoppings
lattice_staggering = -0.75; %SSH-type staggering
raman = 0.3; %Strength of the Raman (flux-carrying) hoppings
RF = 0; %Zero-phase internal transition strength
p = 1; 
q = 4; %Flux per unit cell is 2*pi * (p/q)
num_ks = 100;
k_vals = 2*pi*((1:num_ks) - (num_ks/2))/num_ks;
times = 0:0.1:30;
site1 = 1;
site2 = 32;
current_site = sites/2;
cell_size = 4;
hopping_range = 1;
trap = 0.001; %Trap frequency
quench_strength = 1;
%*********************************************

if mod(site1,q) ~= 1 || mod(site2,q) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

if mod(current_site,q) ~= 0
    error('Current must be measured between unit cells');
end

assert(mod(q,2) == 0,...
    'Unit cell must have even number of sites for SSH staggering');

cells = sites/q;

hoppings_clean = ones(1,2*q)*t;
hoppings_staggered = repmat([t+lattice_staggering/2,t+lattice_staggering/2,...
    t-lattice_staggering/2,t-lattice_staggering/2],1,q/2);
hoppings_imbalanced = repmat([t+laser_detuning/2,t-laser_detuning/2],1,q);
hoppings_both = hoppings_staggered + hoppings_imbalanced - hoppings_clean;

ins1 = TopologicalInsulator_Ladder(sites,hoppings_both,RF,raman,q,p,trap);
ins2 = TopologicalInsulator_Ladder(sites,hoppings_staggered,RF,raman*(1 - quench_strength),q,p,trap);

ins3 = TopologicalInsulator_Ladder(sites,hoppings_staggered,RF,raman,q,p,trap);
ins4 = TopologicalInsulator_Ladder(sites,hoppings_staggered,RF,raman*(1 - quench_strength),q,p,trap);

%ins2 = TopologicalInsulator_C(t,mu,del_2,sites,open);
gap = min(abs(ins1.spectrum));
disp(["Gap = ", num2str(gap)]);
chi = TopologicalInsulator_Ladder.test_symmetries(ins1.hamiltonian);
if abs(chi) > 1.e-3
    disp(["Chial symmetry broken ", num2str(chi)]);
else
    disp("Chial symmetry preserved");
end

init_spinors_1 = ins1.BL_ground_state_spinors(k_vals);
init_spinors_3 = ins3.BL_ground_state_spinors(k_vals);
invar = TopologicalInsulator.BL_wilson_loops(init_spinors_1);
disp(["CS invariant = ", num2str(invar)]);

time_invars_1 = ins2.BL_topological_invariant(init_spinors_1,times,k_vals);
time_invars_3 = ins4.BL_topological_invariant(init_spinors_3,times,k_vals);

charge_op = diag([ones(1,2*current_site),zeros(1,2*(sites-current_site))]);

init_mat_1 = ins1.half_filled_correlation_matrix(-0.1);
init_mat_3 = ins3.half_filled_correlation_matrix(-0.1);
init_charge_1 = real(trace(charge_op * init_mat_1.')); %init_charge_1 = floor(init_charge_1);
init_charge_3 = real(trace(charge_op * init_mat_3.'));
charges_2 = zeros(1,numel(times));
charges_4 = zeros(1,numel(times));

for t_index = 1:numel(times)
    t = times(t_index);
    corrmat_t_2 = ins2.time_evolve_correlation_matrix(init_mat_1,t);
    corrmat_t_4 = ins4.time_evolve_correlation_matrix(init_mat_3,t);
    charges_2(1,t_index) = real(trace(charge_op * corrmat_t_2.')) - init_charge_1;
    charges_4(1,t_index) = real(trace(charge_op * corrmat_t_4.')) - init_charge_3;
end

%% Plotting

time_invars_1 = mod(time_invars_1+0.5,1)-0.5;
charges_2 = mod(charges_2+0.5,1)-0.5;
time_invars_3 = mod(time_invars_3+0.5,1)-0.5;
charges_4 = mod(charges_4+0.5,1)-0.5;

colblue = [0,0.4470,0.7410];
colred = [0.8,0.2,0.05];
lw = 0.75;

jump_inds_1 = [1];
jump_inds_2 = [1];
jump_inds_3 = [1];
jump_inds_4 = [1];
for j = 1:(numel(time_invars_1)-1)
    if abs(time_invars_1(j+1) - time_invars_1(j)) > 0.8
        jump_inds_1 = [jump_inds_1, j,j+1];
    end
    if abs(charges_2(j+1) - charges_2(j)) > 0.8
        jump_inds_2 = [jump_inds_2, j,j+1];
    end
    if abs(time_invars_3(j+1) - time_invars_3(j)) > 0.8
        jump_inds_3 = [jump_inds_3, j,j+1];
    end
    if abs(charges_4(j+1) - charges_4(j)) > 0.8
        jump_inds_4 = [jump_inds_4, j,j+1];
    end
end
jump_inds_1 = [jump_inds_1,numel(time_invars_1)];
jump_inds_2 = [jump_inds_2,numel(time_invars_1)];
jump_inds_3 = [jump_inds_3,numel(time_invars_1)];
jump_inds_4 = [jump_inds_4,numel(time_invars_1)];

figure_handles{end+1} = figure('Name','CS invariant and current');
hold on;
for k = 1:(numel(jump_inds_1)/2)
    st = jump_inds_1(2*k-1); ed = jump_inds_1(2*k);
    h1 = plot(times(st:ed),time_invars_1(st:ed),'LineWidth',lw,'Color',colblue,'DisplayName','$k$-space calculation $(J_1 \neq J_2)$');
end
for k = 1:(numel(jump_inds_2)/2)
    st = jump_inds_2(2*k-1); ed = jump_inds_2(2*k);
    h2 = plot(times(st:ed),-charges_2(st:ed),'--','LineWidth',lw,'Color',colblue,'DisplayName','Experimental simulation $(J_1 \neq J_2)$');
end
for k = 1:(numel(jump_inds_3)/2)
    st = jump_inds_3(2*k-1); ed = jump_inds_3(2*k);
    h3 = plot(times(st:ed),time_invars_3(st:ed),'LineWidth',lw,'Color',colred,'DisplayName','$k$-space calculation $(J_1 = J_2)$');
end
for k = 1:(numel(jump_inds_4)/2)
    st = jump_inds_4(2*k-1); ed = jump_inds_4(2*k);
    h4 = plot(times(st:ed),-charges_4(st:ed),'--','LineWidth',lw,'Color',colred,'DisplayName','Experimental simulation $(J_1 = J_2)$');
end
hold off;
xlabel('Time $t$','interpreter','latex');
ylabel('CS$_1(t)$ mod 1','interpreter','latex');
le = legend([h1,h2,h3,h4]);
le.Interpreter = 'latex';
le.Location = 'NorthWest';
%ylim([-0.05,0.2]);

figure_handles{end+1} = figure('Name','Discrepancy');
plot(times,abs(mod(time_invars_1+charges_2+0.5,1)-0.5));
hold on;
plot(times,times.^2 * (trap * t.^2 / (1)) ); %.*(mod(time_invars+0.5,1)-0.5)
hold off;
set(gca,'Yscale','log');
set(gca,'Xscale','log');