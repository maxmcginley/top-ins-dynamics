classdef MajoranaMemory_Params
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        timestep(1,1) double
        num_steps(1,1) uint32
        steps_per_measure(1,1) uint32
        max_exp(1,1) uint32
        insulator_init cell = {}
        hamiltonians cell = {}
        spectra cell = {}
        reals(1,1) uint32
        num_insulators(1,1) uint32 = 0
        majorana_limit(1,1) double = 0.05
        subspace_plus(2,2) double = 0.5*ones(2,2)
        subspace_minus(2,2) double = 0.5*[[1,-1];[-1,1]]
        spec_freq_widths(1,:) double
        adiabatic_rate double
        spec_amps(1,:) double
        double_noise = false
        system_params
        run_parallel = false
        use_cutoff = false
        cutoffs
        rand_phase = true
    end
    
    methods
        function obj = MajoranaMemory_Params()
        end
        
        function num = num_vals(obj)
            num = obj.num_steps/obj.steps_per_measure;
        end
        
        function times = data_times(obj)
            times = double(0:obj.num_vals()).*(obj.timestep*double(obj.steps_per_measure));
        end
        
        function times = full_times(obj)
            times = double(0:obj.num_steps).*obj.timestep;
        end
    end
    
    methods (Static)
        function outputs = run_jobs(params_in,detailed_save)
            assert(iscell(params_in),'Parameter values to be given as cell');
            num_jobs = numel(params_in);
            
            outputs =  cell(size(params_in));
            
            for j = 1:num_jobs
                maj_params = params_in{j};
                [final_state_minus,final_state_plus] = calculate_majorana_evolution (maj_params);
                
                out_str = struct('maj_params',maj_params);
                
                if detailed_save
                    out_str.final_state_minus = final_state_minus;
                    out_str.final_state_plus = final_state_plus;
                else
                    fids = compute_majorana_fidelities(final_state_minus,final_state_plus);
                    out_str.fidelities = fids;
                end
                outputs{j} = out_str;
            end
        end
    end
end

