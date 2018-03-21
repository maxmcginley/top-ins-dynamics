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
sites = 20;
open = false;
hopping = 1;
nn = 0.3;
phase_1 = 0.4;
phase_2 = -0.1;
sublattice = 0.02;
times = 0:0.1:2;
site1 = 1;
site2 = 40;
cell_size = 4;
hopping_range = 1;
num_ks = 100;
frac_k_window = 1; %Fraction of BZ around zero
%*********************************************

if mod(site1,2) ~= 1 || mod(site2,2) ~= 0
    error('Entanglement cut must not be within a unit cell');
end

k_vals = ([0:(num_ks)] - num_ks/2)*2*pi*frac_k_window/num_ks;

spec = zeros(sites*4,numel(k_vals));

ent_spec = zeros(1+site2-site1,numel(k_vals),numel(times));
ent_spec_aux = zeros(1+site2-site1,numel(k_vals));

gaps_1 = zeros(1,numel(k_vals));
gaps_2 = zeros(1,numel(k_vals));

hamilts = cell(1,numel(k_vals));

for k_index = 1:numel(k_vals)
    ins_k = TopologicalInsulator_Hal(hopping,nn,phase_1,sublattice,k_vals(k_index),sites,open);
    ins_k_2 = TopologicalInsulator_Hal(hopping,nn,phase_2,sublattice,k_vals(k_index),sites,open);
    hamilts{k_index} = ins_k.hamiltonian;
    spec(:,k_index) = ins_k.spectrum;
    
    gaps_1(1,k_index) = min(ins_k.spectrum);
    gaps_2(1,k_index) = min(ins_k_2.spectrum);
    
    corrmat_init = ins_k.half_filled_correlation_matrix(-0.8);
    
    corrmat_aux = ins_k_2.half_filled_correlation_matrix(-0.8);
    for t_index = 1:numel(times)
        time = times(t_index);
        corrmat_t = ins_k_2.time_evolve_correlation_matrix(corrmat_init,time);
        ent_spec(:,k_index,t_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(corrmat_t,site1,site2)));
    end
    ent_spec_aux(:,k_index) = sort(real(TopologicalInsulator.entanglement_spectrum_from_correlation_matrix(corrmat_aux,site1,site2)));
end

%% Gaps

fprintf('Initial gap is %f',min(abs(gaps_1)));
fprintf('Final gap is %f',min(abs(gaps_2)));

%% Plotting

%cont_spec = extract_continuous_spectrum(ent_spec);

figure_handles{end+1} = figure('Name','Kane Mele Entanglement Spectrum');
hold on;
for t_index = 1:numel(times)
    plot(k_vals/(2*pi),cont_spec(:,t_index).','Color',[sqrt(t_index/numel(times)),0,0]);
    %ylim([-0.3,0.3]);
end
hold off;

figure_handles{end+1} = figure('Name','Final Hamiltonian spetrum');
plot(k_vals/(2*pi),spec);
ylim([-4,4]);

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
