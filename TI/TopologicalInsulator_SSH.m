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
            obj@TopologicalInsulator(ham,2);
            obj.hoppingA = hoppingA;
            obj.hoppingB = hoppingB;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
            ham_k = [[0, obj.hoppingA + (obj.hoppingB')*exp(1i*k)];[0,0]];
            ham_k = ham_k + ham_k';
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
        
        function fig_handle = plot_current_entanglement(times,top_invars,deriv_top_invars,r_times,charges,r_currs,entanglements)
            pt_to_inch = 72;
            width = 246; height = 140;
            fig_handle = figure('Name', 'Entanglement spectrum',...
                'Units','inches','Position',[4,3,width/pt_to_inch,height/pt_to_inch]);
            lw = 0.75;
            
            hold on;
            yyaxis left;
            plot(times,top_invars-1,'LineWidth',lw);
            plot(r_times,charges,'--','LineWidth',lw);
            ylim([-0.8,0]);
            ylabel('CS$_1(t) - $ CS$_1(0)$; $Q(t)$','interpreter','latex');
            yyaxis right;
            plot(times(2:(end-1)),deriv_top_invars,'LineWidth',lw);
            plot(r_times,r_currs,'-.','LineWidth',lw);
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

