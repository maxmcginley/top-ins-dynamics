classdef TopologicalInsulator_DoubleKitaev < TopologicalInsulator
    %2 band class BDI Hamiltonian with longer range hoppings
    %   Detailed explanation goes here
    
    properties
        mu
        del
        % t = 1
    end
    
    methods
        function obj = TopologicalInsulator_DoubleKitaev(mu,del,sites,open)
            ham = TopologicalInsulator_DoubleKitaev.DoubleKitaev_hamiltonian(mu,del,sites,open);
            obj@TopologicalInsulator(ham,4);
            obj.mu = mu;
            obj.del = del;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
           error('Not implemented');
        end
        
        %rho_q is the density matrix of the 2-level qubit system
        function corrmat_edge = coherent_majorana_qubit(obj,majorana_limit,rho_q)
            assert(all(size(rho_q) == [2,2]),'qubit must be provided as 2x2 density matrix');
            maj_indices = find(abs(obj.spectrum) < majorana_limit);
            assert(numel(maj_indices) == 4,'Should be 4 Majoranas present');
            
            subspace = kron(rho_q,eye(2));
            
            diag_mat = diag(double(obj.spectrum < 0));
            %diag_mat = eye(numel(obj.spectrum))*0.5;
            diag_mat(maj_indices,maj_indices) = subspace;
            corrmat_edge = obj.orbitals' * diag_mat * obj.orbitals;
        end
    end
    
    methods (Static)
        function ham_full = DoubleKitaev_hamiltonian(mu,del, sites,open)
            cells = sites/4;
            
            diags = [-mu,mu]*0.5;
            hop1 = del*[0,-1];
            hop2 = [-1,1];
            hop3 = del*[1,0];
            
            ham0 = TopologicalInsulator.off_diagonal_matrix(0,diags,cells,open);
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,hop1,cells,open);
            ham2 = TopologicalInsulator.off_diagonal_matrix(2,hop2,cells,open);
            ham3 = TopologicalInsulator.off_diagonal_matrix(3,hop3,cells,open);
            
            ham = ham0 + ham1 + ham2 + ham3;
            ham = ham + ham';
            
            ham_full = [[ham, zeros(size(ham))];[zeros(size(ham)), ham]];
        end
        
        
        
        function U = dirac_to_majorana_matrix(sites)
            if mod(sites,2) ~= 0
                error('matrix not defined on TIs with odd sites');
            end
            U = kron(eye(sites/2),[[1,1];[-1i,1i]]);
        end

        
        %*************SYMMETRIES********************
        
        function C = nambu_operator(index,sites)
            if mod(sites,2) ~= 0
                error('Nambu operator not defined on TIs with odd sites');
            end
            s = cell(3,1);
            s{1} = [[0 1];[1 0]];
            s{2} = [[0 -1i];[1i 0]];
            s{3} = [[1 0];[0 -1]];
            s{4} = eye(2);
            C =kron(eye(sites/2),s{index});
        end
%         
        function phs = test_phs(projector)
            sites = size(projector,1);
            phsop = TopologicalInsulator_DoubleKitaev.nambu_operator(1,sites);
            mat_phs = phsop *  conj(projector) * phsop';

            phs = sum(sum(abs(projector + mat_phs)));
            
        end
    end
end

