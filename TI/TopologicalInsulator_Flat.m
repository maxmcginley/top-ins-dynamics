classdef TopologicalInsulator_Flat < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        intracells
        intercells
    end
    
    methods
        function obj = TopologicalInsulator_Flat(en_scale,winding,sites,open)
            assert(ismember(winding,[-1,0,1]),...
                'Flat band winding can be only 0,1, or -1.');
            [intracells,intercells] = TopologicalInsulator_Flat.hopping_values(...
                en_scale,winding);
            obj@TopologicalInsulator(TopologicalInsulator_Flat.Flat_hamiltonian(...
                intracells,intercells,sites,open));
            obj.intercells = intercells;
            obj.intracells = intracells;
        end
        
        function ham_k = BL_k_hamiltonian(obj,k)
            [t,r,s,u] = deal(obj.intracells);
            [a,b,c,d] = deal(obj.intercells);
            ham_k = [[0, t + a*exp(-1i*k), 0, r + b*exp(1i*k)];...
                     [0, 0, s + c*exp(1i*k), 0];...
                     [0, 0, 0, u + d*exp(1i*k)];...
                     [0,0,0,0]];
            ham_k = ham_k + ham_k';
        end
        
        function [nus,js] = BL_topological_invariant(obj,init_spinors,times,k_vals)
            assert(iscell(init_spinors),'Initial spinors must be provided as a cell');
            assert(numel(init_spinors) == numel(k_vals),...
                'Number of k values must match number of spinors');
            bloch_vectors = zeros(4,numel(k_vals),numel(times));
            for j = 1:numel(k_vals)
                k = k_vals(j);
                hamilt_k = obj.BL_k_hamiltonian(k);
                for t_index = 1:numel(times)
                    time = times(t_index);
                    evol = expm(-1i*hamilt_k*time);
                    bloch_vectors(:,j,t_index) = evol * init_spinors{j};
                end
            end
            
            nus = zeros(1,numel(times));
            js = zeros(1,numel(times));
            
            for t_index = 1:numel(times)
                nus(1,t_index) = TopologicalInsulator_SSH.BL_wilson_loops(...
                    bloch_vectors(:,:,t_index));
                ck = zeros(1,numel(k_vals));
                for j = 1:numel(k_vals)
                    k = k_vals(j);
                    current_k = [[0,1i*(obj.hoppingB)*exp(1i*k)];[0,0]];
                    current_k = current_k + current_k';
                    sp = bloch_vectors(:,j,t_index);
                    ck(j) = sp' * current_k * sp;
                end
                js(1,t_index) = sum(ck)/numel(k_vals);
            end
        end
    end
    
    methods (Static)
        function ham = Flat_hamiltonian(intracells,intercells, sites,open)
            assert(mod(sites,4) == 0,'Must have 4n sites');
            cells = sites/4;
            [t,r,s,u] = deal(intracells);
            [a,b,c,d] = deal(intercells);
            
            hop1 = TopologicalInsulator.off_diagonal_matrix(1,[t',s',u',b],cells,open);
            hop3 = TopologicalInsulator.off_diagonal_matrix(3,[r',0,c,d],cells,open);
            hop5 = TopologicalInsulator.off_diagonal_matrix(5,[a',0,0,0],cells,open);
            
            ham = hop1 + hop3 + hop1' + hop3' + hop5 + hop5';
        end
        
        function [intracells,intercells] = hopping_values(en_scale,winding)
            switch winding
                case 1
                    % t, r, s, u
                    intracells = [1,1i,1i,-1].*en_scale;
                    % a, b, c, d
                    intercells = [1i,-1,1,-1i].*en_scale;
                case -1
                    % t, r, s, u
                    intracells = [1,-1i,-1i,-1].*en_scale;
                    % a, b, c, d
                    intercells = [-1i,-1,1,1i].*en_scale;
                case 0
                    error('Not currently supported');
                otherwise
                    error('Invalid winding');
            end
        end
    end
end

