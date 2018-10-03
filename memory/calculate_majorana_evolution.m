function [final_state_minus,final_state_plus] = calculate_majorana_evolution (maj_params)

state_minus = cell(1,maj_params.num_insulators);
state_plus = cell(1,maj_params.num_insulators);
%final_state_minus = cell(1,maj_params.num_insulators);
%final_state_plus = cell(1,maj_params.num_insulators);

data_times = maj_params.data_times();
full_times = maj_params.full_times();

final_state_minus_rs = zeros([maj_params.reals,maj_params.num_insulators,size(state_minus{1}),numel(data_times)]);
final_state_plus_rs = zeros([maj_params.reals,maj_params.num_insulators,size(state_minus{1}),numel(data_times)]);

for i = 1:maj_params.num_insulators
    display(maj_params.insulator_init{i}.spectrum)
    
    state_minus{i} = maj_params.insulator_init{i}.coherent_majorana_qubit(maj_params.majorana_limit,maj_params.subspace_plus);
    state_plus{i} = maj_params.insulator_init{i}.coherent_majorana_qubit(maj_params.majorana_limit,maj_params.subspace_minus);
end

spectra = TimeEvolution_Noise.generate_poissonians(maj_params.spec_freq_widths,maj_params.spec_amps);
    
tevol = cell(1,maj_params.num_insulators);
    
for i = 1:maj_params.num_insulators
    %final_state_minus(i,:,:,:) = zeros([size(state_minus{i}),numel(data_times)]);
    %final_state_plus(i,:,:,:) = zeros([size(state_plus{i}),numel(data_times)]);
    
    tevol_tmp = TimeEvolution_Noise(maj_params.timestep,maj_params.num_steps,maj_params.max_exp,...
        maj_params.hamiltonians(i,:),spectra,maj_params.reals,false,maj_params.insulator_init{i}.hamiltonian);
    tevol{i} = tevol_tmp.allocate_phases(full_times);
end

prog_handle = waitbar(0,'Time evolving...');  

%% Test symmetries

% [init_trs,init_phs,init_chi] = TopologicalInsulator_DIII.test_symmetries(eye(size(state_plus{1})) - 2*state_plus{1});
% %[init_trs,init_phs,init_chi] = TopologicalInsulator_DIII.test_symmetries(ins_2{1}.hamiltonian);
% 
% [init_trs_BDI,init_phs_BDI,init_chi_BDI] = TopologicalInsulator_ChiralMaj.test_symmetries(eye(size(state_plus{3})) - 2*state_plus{3});
% %[init_trs_BDI,init_phs_BDI,init_chi_BDI] = TopologicalInsulator_ChiralMaj.test_symmetries(ins_1{3}.hamiltonian);
% 
% 
% double_phs = TopologicalInsulator_DoubleKitaev.test_phs(ins_2{2}.hamiltonian);
% double_phs = TopologicalInsulator_DoubleKitaev.test_phs(eye(size(state_plus{2})) - 2*state_plus{2});
% 
% crit = 1.e-6;
% 
% fprintf('DIII: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
% fprintf('D: PHS = %d \n',abs(double_phs) < crit);
% fprintf('BDI: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs_BDI) < crit,abs(init_trs_BDI) < crit,abs(init_chi_BDI) < crit);

%% Time evolve
    reals = double(maj_params.reals);
    
    if maj_params.run_parallel
        num_cores = Inf;
    else
        num_cores = 0;
    end

    parfor (dis_index = 1:maj_params.reals,num_cores)
        
%         if ~maj_params.run_parallel
%             waitbar((dis_index-1)/reals,prog_handle);
%         end
        
        
        
        for i = 1:maj_params.num_insulators
            final_state_minus_real = state_minus{i};
            final_state_plus_real = state_plus{i};
            
%             final_state_minus(dis_index,i,:,:,1) = (final_state_minus_real);
%             final_state_plus(dis_index,i,:,:,1) = (final_state_plus_real);
        
            for t_index = 1:(maj_params.num_vals()+1)

                if t_index > 1
                final_state_minus_real = tevol{i}.evolve(...
                    final_state_minus_real{i},(t_index)*maj_params.steps_per_measure,...
                    maj_params.steps_per_measure,dis_index);
                final_state_plus_real = tevol{i}.evolve(...
                    final_state_plus_real{i},(t_index)*maj_params.steps_per_measure,...
                    maj_params.steps_per_measure,dis_index);
                end
                
                final_state_minus_rs(dis_index,i,:,:,t_index) = (final_state_minus_real);
                final_state_plus_rs(dis_index,i,:,:,t_index) = (final_state_plus_real);
               
            end
        end
        
%         prev_index = max(q_index-1,1);
%         
%         time = TIME_STEP*(q_index - prev_index);
%         
%         final_state_TRS_minus_real =  ...
%             ins_TRS_2.time_evolve_correlation_matrix(final_state_TRS_minus_real,time);
%         final_state_TRS_plus_real =  ...
%             ins_TRS_2.time_evolve_correlation_matrix(final_state_TRS_plus_real,time);
%         final_state_Double_minus_real = ...
%             ins_Double_2.time_evolve_correlation_matrix(final_state_Double_minus_real,time);
%         final_state_Double_plus_real =  ...
%             ins_Double_2.time_evolve_correlation_matrix(final_state_Double_plus_real,time);
%         
%         if mod(q_index - 1,quenches_per_step) == 0
%             t_index = (q_index - 1)/quenches_per_step + 1;
%             
%         end
    end
    
    final_state_minus = reshape(mean(final_state_minus_rs,1),
    

if ishandle(prog_handle)
    close(prog_handle);
end

end

