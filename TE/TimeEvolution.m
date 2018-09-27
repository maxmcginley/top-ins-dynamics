classdef (Abstract) TimeEvolution
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        timestep
        num_steps
        max_exp
    end
    
    methods
        function obj = TimeEvolution(timestep,num_steps,max_exp)
            obj.timestep = timestep;
            obj.num_steps = num_steps;
            obj.max_exp = max_exp;
        end
        
        function state_out = evolve(obj,state_in,step_in,num_steps,varargin)
            if numel(varargin) > 0
                ham_misc = varargin{1};
            else
                ham_misc = [];
            end
            curr_state = state_in;
            ham = obj.calculate_hamiltonian(step_in,ham_misc);
            ev = TimeEvolution.approximate_exponential(1i*ham*obj.timestep,obj.max_exp);
            for s_index = 1:(num_steps/2)
                ham_mid = obj.calculate_hamiltonian(step_in + 2*s_index - 1,ham_misc);
                ham_next = obj.calculate_hamiltonian(step_in + 2*s_index,ham_misc);
                
                ev_mid = TimeEvolution.approximate_exponential(1i*ham_mid*obj.timestep,obj.max_exp);
                ev_next = TimeEvolution.approximate_exponential(1i*ham_next*obj.timestep,obj.max_exp);
                
                inc1 = ev * curr_state * ev' - curr_state;
                inc2 = ev_mid*(curr_state + inc1*(1/2))*ev_mid' - curr_state;
                inc3 = ev_mid*(curr_state + inc2*(1/2))*ev_mid' - curr_state;
                inc4 = ev_next*(curr_state + inc3)*ev_next' - curr_state;
                curr_state = curr_state + (inc1 + inc4 + 2*(inc3 + inc2))/6;
                ev = ev_next;
            end
            state_out = curr_state;
        end
    end
    
    methods(Abstract)
        ham = calculate_hamiltonian(obj,t_step)
    end
    
    methods (Static)
        function inc = approximate_exponential(inc_in,max_exp)
            prev_exp = inc_in;
            inc = eye(size(inc_in));
            for n = 1:max_exp
                inc = inc + prev_exp;
                if n < max_exp
                    prev_exp = prev_exp * inc_in *(1./(n+1));
                end
            end
        end
    end
end

