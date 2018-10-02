classdef TopologicalInsulator_ChiralMaj < TopologicalInsulator
    %2 band class BDI Hamiltonian with longer range hoppings
    %   Detailed explanation goes here
    
    properties
        mu
        del1
        t1
        del2
        t2
    end
    
    methods
        function obj = TopologicalInsulator_ChiralMaj(t1,t2,del1,del2,mu,sites,open,varargin)
            ham = TopologicalInsulator_ChiralMaj.ChiralMaj_hamiltonian(t1,t2,del1,del2,mu,sites,open);
            if numel(varargin) > 0
                diag = varargin{1};
            else
                diag = true;
            end
            obj@TopologicalInsulator(ham,2,diag);
            obj.mu = mu;
            obj.del2 = del2;
            obj.t2 = t2;
            obj.del1 = del1;
            obj.t1 = t1;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
           error('Not implemented');
        end
        
        %rho_q is the density matrix of the 2-level qubit system
        function corrmat_edge = coherent_majorana_qubit(obj,majorana_limit,rho_q)
            assert(all(size(rho_q) == [2,2]),'qubit must be provided as 2x2 density matrix');
            maj_indices = find(abs(obj.spectrum) < majorana_limit);
            assert(numel(maj_indices) == 4,'Should be 4 Majoranas present');
%             maj_wav = obj.orbitals(maj_indices,:);
%             left_indices = sum(abs(maj_wav(:,1:(obj.sites/2))).^2) > 0.5;
%             if left_indices(1) == left_indices(2)
%                 pos_indices = maj_indices([1,3]);
%             else
%                 pos_indices = maj_indices([1,2]);
%             end
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
            
            %sym_orbs = obj.orbitals(maj_indices,:)';
            
            diag_mat = diag(double(obj.spectrum < 0));
            diag_mat(maj_indices,maj_indices) = zeros(4);
            %diag_mat = eye(numel(obj.spectrum))*0.5;
            corrmat_edge = obj.orbitals' * diag_mat * obj.orbitals;
            corrmat_edge = corrmat_edge + sym_orbs * subspace * sym_orbs';
        end
    end
    
    methods (Static)
        function ham = ChiralMaj_hamiltonian(t1,t2,del1,del2,mu,sites,open)
            cells = sites/2;
            
            struc_1 = [0,-1];
            struc_2 = [-1,1];
            struc_3_del1 = [1,0];
            struc_3_del2 = [0,-1];
            struc_4 = [-1,1];
            struc_5 = [1,0];
            
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
            
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,del1,cells,open,2,struc_1);
            ham2 = TopologicalInsulator.off_diagonal_matrix(2,t1,cells,open,2,struc_2);
            ham3_1 = TopologicalInsulator.off_diagonal_matrix(3,del1,cells,open,2,struc_3_del1);
            ham3_2 = TopologicalInsulator.off_diagonal_matrix(3,del2,cells,open,2,struc_3_del2);
            ham4 = TopologicalInsulator.off_diagonal_matrix(4,t2,cells,open,2,struc_4);
            ham5 = TopologicalInsulator.off_diagonal_matrix(5,del2,cells,open,2,struc_5);
            
            ham = ham0 + ham1 + ham2 + ham3_1 + ham3_2 + ham4 + ham5;
            ham = ham + ham';
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
        function [trs,phs,chi] = test_symmetries(projector)
            sites = size(projector,1);
            trsop = TopologicalInsulator_ChiralMaj.nambu_operator(4,sites);
            phsop = TopologicalInsulator_ChiralMaj.nambu_operator(1,sites);
            chiop = TopologicalInsulator_ChiralMaj.nambu_operator(1,sites);
            mat_trs = trsop *  conj(projector) * trsop';
            mat_phs = phsop *  conj(projector) * phsop';
            mat_chi = chiop *  projector * chiop';

            phs = sum(sum(abs(projector + mat_phs)));
            trs = sum(sum(abs(projector - mat_trs)));
            chi = sum(sum(abs(projector + mat_chi)));
            
        end
    end
end

