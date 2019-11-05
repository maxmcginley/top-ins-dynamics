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
%             if numel(varargin) > 0
%                 ham_misc = varargin{1};
%             else
%                 ham_misc = [];
%             end

            %ODE SOLVER
%             sz = size(state_in,2);
%             [~,state_out_ri] = ode45(@(t,y) TimeEvolution.ode_function(t,obj.timestep*y,obj),double([step_in,step_in + num_steps]),...
%                 [real(state_in);imag(state_in)]);
%             state_out = reshape(state_out_ri(end,1:(sz*sz)) + 1i*state_out_ri(end,(sz*sz+1):(2*sz*sz)),[sz,sz]);
            
            %Runge-Kutta
            curr_state = state_in;
            ham = obj.calculate_hamiltonian(step_in,0);
            %ev = TimeEvolution.approximate_exponential(1i*ham*obj.timestep,obj.max_exp);
            for s_index = 1:(num_steps/2)
                ham_mid = obj.calculate_hamiltonian(step_in + 2*s_index - 1,2*s_index - 1);
                ham_next = obj.calculate_hamiltonian(step_in + 2*s_index,2*s_index);
                
                %ev_mid = TimeEvolution.approximate_exponential(1i*ham_mid*obj.timestep,obj.max_exp);
                %ev_next = TimeEvolution.approximate_exponential(1i*ham_next*obj.timestep,obj.max_exp);
                
%                 inc1 = ev * curr_state * ev' - curr_state;
%                 inc2 = ev_mid*(curr_state + inc1*(1/2))*ev_mid' - curr_state;
%                 inc3 = ev_mid*(curr_state + inc2*(1/2))*ev_mid' - curr_state;
%                 inc4 = ev_next*(curr_state + inc3)*ev_next' - curr_state;

                inc1 = -2i*obj.timestep*(ham * curr_state); inc1 = inc1 + inc1';
                inc2 = -2i*obj.timestep*(ham_mid * (curr_state + 0.5*inc1)); inc2 = inc2 + inc2';
                inc3 = -2i*obj.timestep*(ham_mid * (curr_state + 0.5*inc2)); inc3 = inc3 + inc3';
                inc4 = -2i*obj.timestep*(ham_next * (curr_state + inc3)); inc4 = inc4 + inc4';
                curr_state = curr_state + (inc1 + inc4 + 2*(inc3 + inc2))/6;
                ham = ham_next;
                %ev = ev_next;
            end
            unit = obj.calculate_static_unitary(obj.timestep * double(num_steps));
            state_out = unit * curr_state * unit';
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
                    prev_exp = prev_exp * inc_in *(1./double(n+1));
                end
            end
        end
        
        function yp = ode_function(t,y,obj)
            h = obj.calculate_hamiltonian(t,[]);
            sz = sqrt(size(y,1)/2); res = 1:sz; ims = (sz+1):(2*sz);
            y = reshape(y,[2*sz,sz]);
            yp = reshape([real(h) * y(ims,:) + imag(h) * y(res,:) - y(ims,:) * real(h) - y(res,:) * imag(h);...
                real(h) * y(res,:) - imag(h) * y(ims,:) - y(res,:) * real(h) + y(ims,:) * imag(h)], [2*sz*sz,1]);
            %yp = 1i * (h * y - y * h);
        end
    end
end

