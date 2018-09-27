classdef TopologicalInsulator_DoubleKitaev < TopologicalInsulator
    %2 band class BDI Hamiltonian with longer range hoppings
    %   Detailed explanation goes here
    
    properties
        mu
        del
        t
    end
    
    methods
        function obj = TopologicalInsulator_DoubleKitaev(t,mu,del,sites,open,varargin)
            ham = TopologicalInsulator_DoubleKitaev.DoubleKitaev_hamiltonian(t,mu,del,sites,open);
            if numel(varargin) > 0
                diag = varargin{1};
            else
                diag = true;
            end
            obj@TopologicalInsulator(ham,4,diag);
            obj.mu = mu;
            obj.del = del;
            obj.t = t;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
           error('Not implemented');
        end
        
        %rho_q is the density matrix of the 2-level qubit system
        function corrmat_edge = coherent_majorana_qubit(obj,majorana_limit,rho_q)
            assert(all(size(rho_q) == [2,2]),'qubit must be provided as 2x2 density matrix');
            maj_indices = find(abs(obj.spectrum) < majorana_limit);
            assert(numel(maj_indices) == 4,'Should be 4 Majoranas present');
            pos_indices = maj_indices(obj.spectrum(maj_indices) > 0);
            assert(numel(pos_indices) == 2,'Should be 2 Positive-eigenvalue Majoranas');
            
            phsop = TopologicalInsulator_DoubleKitaev.nambu_operator(1,numel(obj.spectrum));
            pos_orbs = obj.orbitals(pos_indices,:)';
            sym_orbs = [pos_orbs, phsop * conj(pos_orbs)]; %a_1 a_2 a_1^\dagger a_2^\dagger
            
            subspace = diag([rho_q(2,2),rho_q(2,2),rho_q(1,1),rho_q(1,1)]);
            subspace(1,4) = -rho_q(1,2);
            subspace(2,3) = rho_q(1,2);
            subspace(3,2) = rho_q(2,1);
            subspace(4,1) = -rho_q(2,1);
            
            sym_orbs = obj.orbitals(maj_indices,:)';
            
            diag_mat = diag(double(obj.spectrum < 0));
            diag_mat(maj_indices,maj_indices) = zeros(4);
            %diag_mat = eye(numel(obj.spectrum))*0.5;
            corrmat_edge = obj.orbitals' * diag_mat * obj.orbitals;
            corrmat_edge = corrmat_edge + sym_orbs * subspace * sym_orbs';
        end
    end
    
    methods (Static)
        function ham = DoubleKitaev_hamiltonian(t,mu,del, sites,open)
            cells = sites/2;
            
            struc_1 = [0,-1];
            struc_2 = [-1,1];
            struc_3 = [1,0];
            
%             hop1 = del*[0,-1];
%             hop2 = [-1,1];
%             hop3 = del*[1,0];
            
            if numel(mu) == 1
                diags = [-mu,mu]*0.5;
                ham0 = TopologicalInsulator.off_diagonal_matrix(0,diags,cells,open);
            else
                mu_structure = [-0.5,0.5];
                ham0 = TopologicalInsulator.off_diagonal_matrix(0,mu,cells,open,2,mu_structure);
            end
            
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,del,cells,open,2,struc_1);
            ham2 = TopologicalInsulator.off_diagonal_matrix(2,t,cells,open,2,struc_2);
            ham3 = TopologicalInsulator.off_diagonal_matrix(3,del,cells,open,2,struc_3);
            
            ham = ham0 + ham1 + ham2 + ham3;
            ham = ham + ham';
            ham((cells+1):end,1:cells) = 0;
            ham(1:cells,(cells+1):end) = 0;
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

