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
        reals(1,1) uint32
        num_insulators(1,1) uint32
        majorana_limit(1,1) double = 0.05
        subspace_plus(2,2) double = 0.5*ones(2,2)
        subspace_minus(2,2) double = 0.5*[[1,-1];[-1,1]]
        spec_freq_widths(1,:) double
        spec_amps(1,:) double
        num_channels(1,1) uint32
        system_params
        run_parallel = false
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
end

