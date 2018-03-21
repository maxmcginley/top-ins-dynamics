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
sites = 24;
open = true;
velocity = 0.25;
mass1 = -0.5;
mass2 = -0.3;
inversion = 0;
chi_breaking = 0.5;
times = [0,1,5];
site1 = 25;
site2 = 72;
cell_size = 4;
hopping_range = 1;
num_ks = 25;
frac_k_window = 0.2; %Fraction of BZ around zero
%*********************************************

if mod(site1,cell_size) ~= 1 || mod(site2,cell_size) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

k_vals = ([0:(num_ks)] - num_ks/2)*2*pi*frac_k_window/num_ks;

spec = zeros(numel(k_vals),numel(k_vals),sites*4);

ent_spec = zeros(numel(k_vals),numel(k_vals),1+site2-site1,numel(times));
% ent_spec_aux = zeros(1+site2-site1,numel(k_vals));

gaps_1 = zeros(numel(k_vals),numel(k_vals));
gaps_2 = zeros(numel(k_vals),numel(k_vals));

hamilts = cell(numel(k_vals),numel(k_vals));

for kx_index = 1:numel(k_vals)
    for ky_index = 1:numel(k_vals)
        ins_k = TopologicalInsulator_CII(velocity,mass1,inversion,chi_breaking,k_vals(kx_index),k_vals(ky_index),sites,open);
        ins_k_2 = TopologicalInsulator_CII(velocity,mass2,inversion,chi_breaking,k_vals(kx_index),k_vals(ky_index),sites,open);
        hamilts{kx_index,ky_index} = ins_k.hamiltonian;
        spec(kx_index,ky_index,:) = ins_k.spectrum;

        gaps_1(kx_index,ky_index) = min(abs(ins_k.spectrum));
        gaps_2(kx_index,ky_index) = min(abs(ins_k_2.spectrum));

        corrmat_init = ins_k.half_filled_correlation_matrix(0);
        
        for t_index = 1:numel(times)
            time = times(t_index);
            corrmat_t = ins_k_2.time_evolve_correlation_matrix(corrmat_init,time);
            ent_spec(kx_index,ky_index,:,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(corrmat_t,site1,site2)));
        end
    end
end

%% Test Symmetries

init_chi = TopologicalInsulator_CII.test_chiral_symmetry(hamilts{1,1});

crit = 1.e-6;

fprintf('Hamiltonian: CHI = %d \n',abs(init_chi) < crit);

%% Gaps

fprintf('Initial gap is %f',min(min(abs(gaps_1))));
fprintf('Final gap is %f',min(min(abs(gaps_2))));

%% Plotting

%cont_spec = extract_continuous_spectrum(ent_spec);



figure_handles{end+1} = figure('Name','Class CII Entanglement Spectrum before evolution');
hold on;
for j= 1:(1+site2-site1)
    surf(k_vals/(2*pi),k_vals/(2*pi),ent_spec(:,:,j,1));
end
hold off;
view(17,22);
zlim([-0.8,0.8]);

figure_handles{end+1} = figure('Name','Class CII Entanglement Spectrum after evolution');
hold on;
for j= 1:(1+site2-site1)
    surf(k_vals/(2*pi),k_vals/(2*pi),ent_spec(:,:,j,2));
end
hold off;
view(17,22);
zlim([-0.8,0.8]);

figure_handles{end+1} = figure('Name','Class CII Entanglement Spectrum after long evolution');
hold on;
for j= 1:(1+site2-site1)
    surf(k_vals/(2*pi),k_vals/(2*pi),ent_spec(:,:,j,3));
end
hold off;
view(17,22);
zlim([-0.8,0.8]);

%% Functions

function ents_continuous = extract_continuous_spectrum(ents_in)
    spec_prev = [];
    ents_continuous = zeros(size(ents_in,2),size(ents_in,3));
    for t_index = 1:size(ents_in,3)
        all_specs = ents_in(:,:,t_index);
        if isempty(spec_prev)
            [~,ind] = min(min(abs(all_specs),[],2),[],1);
            disp(ind);
            ents_continuous(:,t_index) = all_specs(ind,:);
        else
            sum_diffs = NaN(1,size(all_specs,1));
            for s_index = 1:size(all_specs,1)
                sum_diffs(1,s_index) = sum(abs(all_specs(s_index,:) - spec_prev).^2);
            end
            [~,ind] = min(sum_diffs);
            ents_continuous(:,t_index) = all_specs(ind,:);
        end
        spec_prev = ents_continuous(:,t_index).';
    end
end
