classdef TopologicalInsulator_SSH < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hoppingA
        hoppingB
    end
    
    methods
        function obj = TopologicalInsulator_SSH(hoppingA,hoppingB,sites,open)
            ham = TopologicalInsulator_SSH.SSH_hamiltonian(hoppingA,hoppingB,sites,open);
            obj@TopologicalInsulator(ham);
            obj.hoppingA = hoppingA;
            obj.hoppingB = hoppingB;
        end
        
        function [nus,js] = BL_topological_invariant(obj,init_spinors,times,k_vals)
            assert(iscell(init_spinors),'Initial spinors must be provided as a cell');
            assert(numel(init_spinors) == numel(k_vals),...
                'Number of k values must match number of spinors');
            bloch_vectors = zeros(2,numel(k_vals),numel(times));
            for j = 1:numel(k_vals)
                k = k_vals(j);
                hamilt_k = [[0,obj.hoppingA + (obj.hoppingB')*exp(1i*k)];[0,0]];
                hamilt_k = hamilt_k + hamilt_k';
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
        function ham = SSH_hamiltonian(hopping1, hopping2, sites,open)
            if mod(sites,2) == 0
                ham = diag([repmat([hopping1, hopping2],1,(sites-2)/2), hopping1], 1);
                if ~open
                    ham(sites,1) = hopping2;
                end
                
                ham = ham + ham';
            else
                if ~open
                    error('Must have even number of sites for periodic system');
                end
                ham = diag(repmat([hopping1, hopping2],1,(sites-1)/2), 1);
                ham = ham + ham';
            end
        end
        
        function sps = BL_constant_spinor(spinor,k_vals)
            sps = cell(1,numel(k_vals));
            sps(:) = {spinor};
        end
        
        function wilson_loop = BL_wilson_loops(spins)
            els = zeros(1,size(spins,2));
            for j = 1:size(spins,2)
                if j ~= size(spins,2)
                    next = j+1;
                else
                    next = 1;
                end
                els(j) = spins(:,j)' * spins(:,next);
            end
            wilson_loop = mod(sum(angle(els))/(2*pi),1);
        end
        
        %*************SYMMETRIES********************
        
        function C = nambu_operator(index,sites)
            if mod(sites,2) ~= 0
                error('Nambu operator not defined on TIs with odd sites');
            end
            switch index
                case 1
                    tau = [[0 1];[1 0]];
                case 2
                    tau = [[0 -1i];[1i 0]];
                case 3
                    tau = [[1 0];[0 -1]];
                otherwise
                    error('Invalid Pauli index');
            end
            C = kron(eye(sites/2),tau);
        end
        
        function [trs,phs,chi] = test_symmetries(corrmat)
            sites = size(corrmat,1);
            mat_phs = TopologicalInsulator_SSH.nambu_operator(1,sites) * ...
                conj(corrmat) * TopologicalInsulator_SSH.nambu_operator(1,sites);
            mat_trs = conj(corrmat);
            mat_chi = TopologicalInsulator_SSH.nambu_operator(3,sites) * ...
                corrmat * TopologicalInsulator_SSH.nambu_operator(3,sites);

            phs = sum(sum(abs((eye(sites) - 2*corrmat) + (eye(sites) - 2*mat_phs))));
            trs = sum(sum(abs((eye(sites) - 2*corrmat) - (eye(sites) - 2*mat_trs))));
            chi = sum(sum(abs((eye(sites) - 2*corrmat) + (eye(sites) - 2*mat_chi))));
            
        end
    end
end

