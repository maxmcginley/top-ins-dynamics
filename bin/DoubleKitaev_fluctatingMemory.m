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

%TCM Computer
%file_directory = fullfile('/','rscratch','mm2025','data','majorana');
%Local
file_directory = '.';

filename = fullfile(file_directory,'majorana_memory_out.mat');

detailed_save = false;

VARY_AMPLITUDE = true;
TESTING = false;

%******************INPUT DATA*******************
if TESTING
    run(fullfile(file_directory,'memory_input_testing.m'));
else
    run(fullfile(file_directory,'memory_input.m'));
end

RUN = false;
SAVE = true;
%*********************************************

if RUN
    
    
    while true
        reply = input('Running new simulation... Previous data saved?','s');
        if strcmp(reply,'q') || strcmp(reply,'quit')
            disp('Exiting...');
            return
        elseif strcmp(reply,'yes')
            break
        else
            disp('Unrecognised input\n');
        end
    end
    
    profile off;
    profile on;
    
    if VARY_AMPLITUDE
        params_in = cell(size(amplitudes,1),1);
        for amp_index = 1:size(amplitudes,1)
            maj_copy = maj_params;
            maj_copy.spec_amps = amplitudes(amp_index,:);
            maj_copy = setup_insulators(maj_copy);
            params_in{amp_index,1} = maj_copy;
        end
    else
        maj_params = setup_insulators(maj_params);
        params_in = {maj_params};
    end

    outputs = MajoranaMemory_Params.run_jobs(params_in,detailed_save);
    
    profile off;
    profile viewer;
end
%% Saving data

if RUN && SAVE
    save(filename,'outputs');
end

%% Plotting

if ~RUN
outputs_struct = load(filename);
outputs = outputs_struct.outputs;

if numel(outputs) == 1
    output = outputs{1};
    if isfield(output,'final_state_minus') && isfield(output,'final_state_plus')
        fidelities = compute_majorana_fidelities(output.final_state_minus,output.final_state_plus);
    elseif isfield(output,'fidelities')
        fidelities = output.fidelities;
    else
        error('Invalid output format')
    end

    figure_handles{end+1} = majorana_memory_plot(fidelities,output.maj_params.data_times());
else
    figure_handles{end+1} = majorana_amplitude_plot(outputs);
end

end
