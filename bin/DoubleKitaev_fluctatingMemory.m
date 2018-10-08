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
addpath(fullfile(pwd,'..','TE'));
addpath(fullfile(pwd,'..','memory'));

filename = 'majorana_memory_out.mat';

%******************INPUT DATA*******************
maj_params = MajoranaMemory_Params();

maj_params.timestep = 0.1;
maj_params.num_steps = 8000;
maj_params.steps_per_measure = 40;

maj_params.max_exp = 3;

maj_params.reals = 50;

maj_params.run_parallel = true;

system_params = struct();

system_params.sites = 16;
system_params.t = 1;
system_params.mu_init = 0.25;
system_params.del_1 = 1;
system_params.del_2 = 0.3;
system_params.alpha = 0.2;
system_params.mu_patch_size = 2;

maj_params.system_params = system_params;

maj_params.spec_freq_widths = [400,400];
maj_params.spec_amps = [0.025,0.0125];

RUN = false;
SAVE = true;
%*********************************************

if RUN
    maj_params = setup_insulators(maj_params);    

    [final_state_minus,final_state_plus] = calculate_majorana_evolution (maj_params);

    fidelities = compute_majorana_fidelities(final_state_minus,final_state_plus);
end
%% Saving data

if RUN && SAVE
    output = struct('final_state_minus',{final_state_minus},'final_state_plus',{final_state_plus},...
        'maj_params',maj_params,'fidelities',{fidelities});
    save(filename,'-struct','output');
end

%% Plotting

if ~RUN
output = load(filename);

figure_handles{end+1} = majorana_memory_plot(output.fidelities,output.maj_params.data_times());

end
