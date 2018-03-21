classdef TopologicalInsulator_DIII < TopologicalInsulator
    %2 band class BDI Hamiltonian with longer range hoppings
    %   Detailed explanation goes here
    
    properties
        mu
        delp
        dels
        alpha
    end
    
    methods
        function obj = TopologicalInsulator_DIII(mu,delp,dels,alpha,sites,open)
            ham = TopologicalInsulator_DIII.DIII_hamiltonian(mu,delp,dels,alpha,sites,open);
            obj@TopologicalInsulator(ham,4);
            obj.mu = mu;
            obj.delp = delp;
            obj.dels = dels;
            obj.alpha = alpha;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
           error('Not implemented');
        end
        
        
    end
    
    methods (Static)
        function ham = DIII_hamiltonian(mu,delp,dels,alpha, sites,open)
            cells = sites/4;
            
            diags = [1-mu,1-mu,mu-1,mu-1];
            hop1 = dels*[0,1,0,0];
            hop2 = delp*[-1,-1,0,0];
            hop3 = [-dels,1i*alpha,0,-1i*alpha];
            hop4 = [-1,-1,1,1]*0.5;
            hop5 = [1i*alpha,0,-1i*alpha,0];
            hop6 = delp*[0,0,1,1];
            
            ham0 = TopologicalInsulator.off_diagonal_matrix(0,diags,cells,open);
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,hop1,cells,open);
            ham2 = TopologicalInsulator.off_diagonal_matrix(2,hop2,cells,open);
            ham3 = TopologicalInsulator.off_diagonal_matrix(3,hop3,cells,open);
            ham4 = TopologicalInsulator.off_diagonal_matrix(4,hop4,cells,open);
            ham5 = TopologicalInsulator.off_diagonal_matrix(5,hop5,cells,open);
            ham6 = TopologicalInsulator.off_diagonal_matrix(6,hop6,cells,open);
            
            ham = ham0 + ham1 + ham2 + ham3 + ham4 + ham5 + ham6;
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
%         
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

