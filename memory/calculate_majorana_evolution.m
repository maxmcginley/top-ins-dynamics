function [final_state_minus,final_state_plus] = calculate_majorana_evolution (maj_params)

state_minus = cell(1,maj_params.num_insulators);
state_plus = cell(1,maj_params.num_insulators);
%final_state_minus = cell(1,maj_params.num_insulators);
%final_state_plus = cell(1,maj_params.num_insulators);

data_times = maj_params.data_times();
full_times = maj_params.full_times();

final_state_minus_rs = cell(maj_params.reals,maj_params.num_insulators);
final_state_plus_rs = cell(maj_params.reals,maj_params.num_insulators);


for i = 1:maj_params.num_insulators
    display(['Min energy is ',num2str(min(abs(maj_params.insulator_init{i}.spectrum)))]);
    
    state_minus{i} = maj_params.insulator_init{i}.coherent_majorana_qubit(maj_params.majorana_limit,maj_params.subspace_plus);
    state_plus{i} = maj_params.insulator_init{i}.coherent_majorana_qubit(maj_params.majorana_limit,maj_params.subspace_minus);
    for dis_index = 1:maj_params.reals
        final_state_minus_rs{dis_index,i} = NaN([size(state_minus{i}),numel(data_times)]);
        final_state_plus_rs{dis_index,i} = NaN([size(state_plus{i}),numel(data_times)]);
    end
    
    
end

test_symmetries(state_plus, maj_params.num_insulators);

%TENTATIVE - explicitly set majornana mode energies to zero
corr_spec = maj_params.insulator_init{1}.spectrum;
corr_spec(abs(corr_spec) < 1.e-6) = 0;
maj_params.insulator_init{1}.spectrum = corr_spec;
maj_params.insulator_init{1}.hamiltonian = maj_params.insulator_init{1}.orbitals' ...
    * diag(maj_params.insulator_init{1}.spectrum) * maj_params.insulator_init{1}.orbitals;
maj_params.insulator_init{1}.hamiltonian = (maj_params.insulator_init{1}.hamiltonian + maj_params.insulator_init{1}.hamiltonian')/2;
%**********************

%% Time evolve
    reals = maj_params.reals;
    
    if maj_params.run_parallel
        num_cores = 3;
    else
        num_cores = 0;
    end
    
    num_ins = maj_params.num_insulators;
    num_vals = maj_params.num_vals();
    steps_per_measure = maj_params.steps_per_measure;
    timestep = maj_params.timestep;
    num_steps = maj_params.num_steps;
    max_exp = maj_params.max_exp;
    hamiltonians = maj_params.hamiltonians;
    ins_init = maj_params.insulator_init;
    spectra = maj_params.spectra;
    
    %**********DEBUG*******
    DEBUG_deterministic = false;
    %**********************

    parfor (dis_index = 1:reals,num_cores)
    %*****DEBUG_MODE***********
    %for dis_index = 1:reals
    %**************************
        
        display(['Running job ',num2str(dis_index)]);
        
        for i = 1:num_ins
            
            tevol_tmp = TimeEvolution_Noise(timestep,num_steps,max_exp,...
            hamiltonians{i},spectra{i},maj_params.rand_phase,ins_init{i}.hamiltonian,true,maj_params.adiabatic_rate,DEBUG_deterministic);
                tevol = tevol_tmp.allocate_phases(full_times);
            
            
                
            final_state_minus_real = tevol.static_evec' * state_minus{i} * tevol.static_evec;
            final_state_plus_real = tevol.static_evec' * state_plus{i} * tevol.static_evec;
            
%             final_state_minus(dis_index,i,:,:,1) = (final_state_minus_real);
%             final_state_plus(dis_index,i,:,:,1) = (final_state_plus_real);
        
            for t_index = 1:(num_vals+1)

                if t_index > 1
                final_state_minus_real = tevol.evolve(...
                    final_state_minus_real,(t_index-2)*steps_per_measure,...
                    steps_per_measure,dis_index);
                final_state_plus_real = tevol.evolve(...
                    final_state_plus_real,(t_index-2)*steps_per_measure,...
                    steps_per_measure,dis_index);
                end
                
                final_state_minus_rs{dis_index,i}(:,:,t_index) = (final_state_minus_real);
                final_state_plus_rs{dis_index,i}(:,:,t_index) = (final_state_plus_real);
               
            end
            
            if DEBUG_deterministic
                stat_ham = tevol.static_ham;
                for k = 1:tevol.num_channels
                    stat_ham = stat_ham + imag(tevol.fourier_amplitudes(k,1))*(tevol.static_evec * tevol.hams{k} * tevol.static_evec');
                end
                time = double(num_steps) * tevol.timestep;
                unit = expm(-1i*stat_ham*time);
                exact = unit * state_minus{i} * unit';
                rk = tevol.static_evec * final_state_minus_real * tevol.static_evec';
                %ph = tevol.calculate_static_phases(time);
                
                display(sum(sum(abs( (exact) - rk))));
            end
            
            es = eig(final_state_minus_real);
            display(['Entropy of final state is ',num2str(-sum(log(es).*es))]);
            
        end
        
    end
    
    final_state_minus = cell(1,num_ins);
    final_state_plus = cell(1,num_ins);
    
    for i=1:num_ins
        final_state_minus{i} = mean(cat(4,final_state_minus_rs{:,i}),4);
        final_state_plus{i} = mean(cat(4,final_state_plus_rs{:,i}),4);
    end
    

% if ishandle(prog_handle)
%     close(prog_handle);
% end

end

