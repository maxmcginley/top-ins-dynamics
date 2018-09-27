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
open = false;

hopping = 1;
so_1 = 0.5;
so_2 = 1.5;
rashba = 0.1;
sublattice = 0.2;

nn = 0.3;
phase_1 = 0.4;
phase_2 = -0.2;
sublattice_hal = 0.1;

times = [0,2];

site1_km = 1;
site2_km = 160;

site1_hal = 1;
site2_hal = 80;

cell_size = 8;
hopping_range = 1;
num_ks = 200;
frac_k_window = 0.2; %Fraction of BZ around zero
%*********************************************

if mod(site1_km,2) ~= 1 || mod(site2_km,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

k_vals = ([0:(num_ks)] - num_ks/2)*2*pi*frac_k_window/num_ks;

spec = zeros(sites*8,numel(k_vals));

ent_spec_hal = zeros(1+site2_hal-site1_hal,numel(k_vals),numel(times));
ent_spec_km = zeros(1+site2_km-site1_km,numel(k_vals),numel(times));

gaps_km = zeros(1,numel(k_vals));
gaps_hal = zeros(1,numel(k_vals));

for k_index = 1:numel(k_vals)
     ins_k_hal = TopologicalInsulator_Hal(hopping,nn,phase_1,sublattice_hal,k_vals(k_index),sites,open);
    ins_k_2_hal = TopologicalInsulator_Hal(hopping,nn,phase_2,sublattice_hal,k_vals(k_index),sites,open);
    
    ins_k_km = TopologicalInsulator_KM(hopping,so_1,rashba,sublattice,k_vals(k_index),sites,open);
    ins_k_2_km = TopologicalInsulator_KM(hopping,so_2,rashba,sublattice,k_vals(k_index),sites,open);
    
    spec(:,k_index) = ins_k_km.spectrum;
    
    gaps_km(1,k_index) = min(ins_k_km.spectrum);
    gaps_hal(1,k_index) = min(ins_k_hal.spectrum);
    
    corrmat_init_hal = ins_k_hal.half_filled_correlation_matrix(-0.8);
    corrmat_init_km = ins_k_km.half_filled_correlation_matrix(-0.8);
    
    for t_index = 1:numel(times)
        time = times(t_index);
        corrmat_t_km = ins_k_2_km.time_evolve_correlation_matrix(corrmat_init_km,time);
        corrmat_t_hal = ins_k_2_hal.time_evolve_correlation_matrix(corrmat_init_hal,time);
        ent_spec_km(:,k_index,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(corrmat_t_km,site1_km,site2_km)));
        ent_spec_hal(:,k_index,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(corrmat_t_hal,site1_hal,site2_hal)));
    end
end

ind_choice = 5;

%TopologicalInsulator_KM.matrix_has_trs(hamilts_km{1,(num_ks/2) + 1 + ind_choice},hamilts_km{1,(num_ks/2) + 1 - ind_choice})

%% Gaps

fprintf('KM gap is %f',min(abs(gaps_km)));
fprintf('Haldane gap is %f',min(abs(gaps_hal)));

%% Plotting

%cont_spec = extract_continuous_spectrum(ent_spec);

%figure_handles{end+1} = figure('Name','Kane Mele Entanglement Spectrum');

figure_handles{end+1} = TopologicalInsulator_KM.entanglement_plot(k_vals,ent_spec_hal(:,:,1),ent_spec_hal(:,:,end),ent_spec_km(:,:,1),ent_spec_km(:,:,end));

% hold on;
% for t_index = 1:numel(times)
%     plot(k_vals/(2*pi),cont_spec(:,t_index).','Color',[sqrt(t_index/numel(times)),0,0]);
%     ylim([-0.3,0.3]);
% end
% hold off;

% figure_handles{end+1} = figure('Name','Final Hamiltonian entanglement spetrum');
% plot(k_vals/(2*pi),ent_spec_aux.');
% ylim([-4,4]);

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
