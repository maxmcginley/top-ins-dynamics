function fig_handle = majorana_amplitude_plot(outputs)

num_outputs = numel(outputs);

decay_consts = NaN(size(outputs));
sat_vals = NaN(size(outputs));
amps = NaN(size(outputs));

for j = 1:num_outputs
    output = outputs{j};
    if isfield(output,'final_state_minus') && isfield(output,'final_state_plus')
        fidelities = compute_majorana_fidelities(output.final_state_minus,output.final_state_plus);
    elseif isfield(output,'fidelities')
        fidelities = output.fidelities;
    else
        error('Invalid output format')
    end
    
    times = output.maj_params.data_times();
    
    y_data = reshape(-log(fidelities{1}/2),numel(fidelities{1}),1);
    reg_data = [ones(numel(fidelities{1}),1), reshape(times,numel(times),1)];
    
    fit_indices = y_data < 10;
    
    regs = reg_data(fit_indices,:) \ y_data(fit_indices);
    decay_consts(j) = regs(2);
    
    sat_vals(j) = 0.5*fidelities{2}(end);
    
    amps(j) = output.maj_params.spec_amps(1);
    
end



width = 246;
height = 150;

pos = [300,200,width,height];

fig_handle = figure('Name','Amplitude scaling','Units','points','Position',pos);
hold on;
plot(amps,decay_consts);
yyaxis right;
plot(amps,sat_vals);
hold off;

set(gca,'XScale','log');
yyaxis left;
set(gca,'YScale','log');

reg_data_decay = [ones(numel(amps),1), reshape(log(amps),numel(amps),1)];
y_data_decay = reshape(log(decay_consts),numel(decay_consts),1);
regs_decay = reg_data_decay \ y_data_decay;

disp(regs_decay(2));

end