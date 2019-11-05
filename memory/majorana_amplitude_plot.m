function fig_handle = majorana_amplitude_plot(outputs)

num_outputs = numel(outputs);

decay_consts = NaN(size(outputs));
sat_vals = NaN(size(outputs));
amps = NaN(size(outputs));

min_x = 1;
max_x = numel(outputs{1}.maj_params.data_times());
%max_x = 7;
min_x = 8;
num_x = max_x - min_x + 1;

fidelities = cell(1,num_outputs);

for j = 1:num_outputs
    output = outputs{j};
    if isfield(output,'final_state_minus') && isfield(output,'final_state_plus')
        fidelities{j} = compute_majorana_fidelities(output.final_state_minus,output.final_state_plus);
    elseif isfield(output,'fidelities')
        fidelities{j} = output.fidelities;
    else
        error('Invalid output format')
    end
    
    times = output.maj_params.data_times();
    times = times(min_x:max_x);
    
    %Logarithmic scale
    y_data = reshape(-log(fidelities{j}{1}(min_x:max_x)/2),num_x,1);
    %Linear scale
    y_data = reshape(-fidelities{j}{1}(min_x:max_x)/2,num_x,1);
    
    reg_data = [ones(num_x,1), reshape(times,num_x,1)];
    
    fit_indices = y_data < 10;
    
    regs = reg_data(fit_indices,:) \ y_data(fit_indices);
    decay_consts(j) = regs(2);
    
    if output.maj_params.num_insulators >= 2
        sat_vals(j) = 0.5*fidelities{j}{2}(end);
    end
    
    amps(j) = output.maj_params.spec_amps(1);
    
end



width = 246;
height = 150;

pos = [300,200,width,height];

figure('Name','Memory decay','Units','points','Position',pos);
hold on;
for j = 1:num_outputs
    plot(times,fidelities{j}{1}(min_x:max_x));
end
hold off;

fig_handle = figure('Name','Amplitude scaling','Units','points','Position',pos);
hold on;
plot(amps,decay_consts);
set(gca,'XScale','log');
set(gca,'YScale','log');
if output.maj_params.num_insulators >= 2
yyaxis right;
plot(amps,sat_vals);
yyaxis left;
end
hold off;

reg_data_decay = [ones(numel(amps),1), reshape(log(amps),numel(amps),1)];
y_data_decay = reshape(log(decay_consts),numel(decay_consts),1);
regs_decay = reg_data_decay \ y_data_decay;

disp(regs_decay(2));

end