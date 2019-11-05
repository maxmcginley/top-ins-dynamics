classdef TopologicalInsulator_DIII < TopologicalInsulator
    %2 band class BDI Hamiltonian with longer range hoppings
    %   Detailed explanation goes here
    
    properties
        mu
        delp
        dels
        alpha
        t
    end
    
    methods
        function obj = TopologicalInsulator_DIII(t,mu,delp,dels,alpha,sites,open,varargin)
            ham = TopologicalInsulator_DIII.DIII_hamiltonian(t,mu,delp,dels,alpha,sites,open);
            if numel(varargin) > 0
                diag = varargin{1};
            else
                diag = true;
            end
            obj@TopologicalInsulator(ham,4,diag);
            obj.mu = mu;
            obj.delp = delp;
            obj.dels = dels;
            obj.alpha = alpha;
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
            
            phsop = TopologicalInsulator_DIII.nambu_operator(1,4,numel(obj.spectrum));
            pos_orbs = obj.orbitals(pos_indices,:)';
            sym_orbs = [pos_orbs, phsop * conj(pos_orbs)]; %a_1 a_2 a_1^\dagger a_2^\dagger 
            
            subspace = diag([rho_q(2,2),rho_q(2,2),rho_q(1,1),rho_q(1,1)]);
            subspace(1,4) = -rho_q(1,2);
            subspace(2,3) = rho_q(1,2);
            subspace(3,2) = rho_q(2,1);
            subspace(4,1) = -rho_q(2,1);
            
            diag_mat = diag(double(obj.spectrum < 0));
            diag_mat(maj_indices,maj_indices) = zeros(4);
            %diag_mat = eye(numel(obj.spectrum))*0.5;
            corrmat_edge = obj.orbitals' * diag_mat * obj.orbitals;
            corrmat_edge = corrmat_edge + sym_orbs * subspace * sym_orbs';
        end
    end
    
    methods (Static)
        function ham = DIII_hamiltonian(t,mu,delp,dels,alpha, sites,open)
            cells = sites/4;
            
            struc_1 = [0,1,0,0];
            struc_2 = [0,0,-1,-1];
            struc_3_dels = [-1,0,0,0];
            struc_3_alpha = [0,1,0,-1];
            struc_4 = [-1,-1,1,1];
            struc_5 = [-1,0,1,0];
            struc_6 = [1,1,0,0];
            
%             hop1 = dels*[0,1,0,0];
%             hop2 = delp*[0,0,-1,-1];
%             hop3 = [-dels,alpha,0,-alpha];
%             hop4 = [-1,-1,1,1];
%             hop5 = [-alpha,0,alpha,0];
%             hop6 = delp*[1,1,0,0];
            
            if numel(mu) == 1
                diags = [-mu,-mu,mu,mu]*0.5;
                ham0 = TopologicalInsulator.off_diagonal_matrix(0,diags,cells,open);
            else
                mu_structure = [-0.5,-0.5,0.5,0.5];
                ham0 = TopologicalInsulator.off_diagonal_matrix(0,mu,cells,open,4,mu_structure);
            end
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,dels,cells,open,4,struc_1);
            ham2 = TopologicalInsulator.off_diagonal_matrix(2,delp,cells,open,4,struc_2);
            ham3_dels = TopologicalInsulator.off_diagonal_matrix(3,dels,cells,open,4,struc_3_dels);
            ham3_alpha = TopologicalInsulator.off_diagonal_matrix(3,alpha,cells,open,4,struc_3_alpha);
            ham4 = TopologicalInsulator.off_diagonal_matrix(4,t,cells,open,4,struc_4);
            ham5 = TopologicalInsulator.off_diagonal_matrix(5,alpha,cells,open,4,struc_5);
            ham6 = TopologicalInsulator.off_diagonal_matrix(6,delp,cells,open,4,struc_6);
            
            ham = ham0 + ham1 + ham2 + ham3_dels + ham3_alpha + ham4 + ham5 + ham6;
            ham = ham + ham';
        end
        
        function fig_handle = plot_current_entanglement(times,top_invars,deriv_top_invars,r_times,charges,r_currs,entanglements)
            pt_to_inch = 72;
            width = 246; height = 170;
            fig_handle = figure('Name', 'Entanglement spectrum',...
                'Units','inches','Position',[4,3,width/pt_to_inch,height/pt_to_inch]);
            
            hold on;
            yyaxis left;
            plot(times,top_invars-1);
            plot(r_times,charges,'s');
            ylim([-0.8,0]);
            ylabel('CS$_1(t) - $ CS$_1(0)$','interpreter','latex');
            yyaxis right;
            plot(times(2:(end-1)),deriv_top_invars,'--');
            plot(r_times,r_currs,'.','MarkerSize',12);
            ylim([-inf,0.3]);
            hold off;
            xlabel('Time $t$','interpreter','latex');
            ylabel('Current $j(t)$','interpreter','latex');
            ax1 = gca;
            
            set(ax1,'TickLabelInterpreter','latex');
            set(ax1,'FontSize',8);
            set(ax1,'FontName','Times');
            
            box on;
        end
        
        function U = dirac_to_majorana_matrix(sites)
            if mod(sites,4) ~= 0
                error('matrix not defined on TIs with odd sites');
            end
            U = kron(eye(sites/4),kron([[1,1];[-1i,1i]],eye(2)));
        end
        
        

        
        %*************SYMMETRIES********************
        
        function C = nambu_operator(index1,index2,sites)
            if mod(sites,2) ~= 0
                error('Nambu operator not defined on TIs with odd sites');
            end
            s = cell(3,1);
            s{1} = [[0 1];[1 0]];
            s{2} = [[0 -1i];[1i 0]];
            s{3} = [[1 0];[0 -1]];
            s{4} = eye(2);
            C = kron(eye(sites/4),kron(s{index1},s{index2}));
        end
%         
        function [trs,phs,chi] = test_symmetries(projector)
            sites = size(projector,1);
            trsop = TopologicalInsulator_DIII.nambu_operator(4,2,sites);
            phsop = TopologicalInsulator_DIII.nambu_operator(1,4,sites);
            chiop = TopologicalInsulator_DIII.nambu_operator(1,2,sites);
            mat_trs = trsop *  conj(projector) * trsop';
            mat_phs = phsop *  conj(projector) * phsop';
            mat_chi = chiop *  projector * chiop';

            phs = sum(sum(abs(projector + mat_phs)))/numel(projector);
            trs = sum(sum(abs(projector - mat_trs)))/numel(projector);
            chi = sum(sum(abs(projector + mat_chi)))/numel(projector);
            
        end
    end
end

