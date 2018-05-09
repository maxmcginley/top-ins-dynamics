classdef TopologicalInsulator_SSH < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hoppingA
        hoppingB
        hopping_A_3
        hopping_B_3
    end
    
    methods
        function obj = TopologicalInsulator_SSH(hoppingA,hoppingB,sites,open,varargin)
            if numel(varargin) > 0
                assert(numel(varargin) >= 2,'Must provide two third-order hoppings');
                hopping_A_3 = varargin{1};
                hopping_B_3 = varargin{2};
            else
                hopping_A_3 = 0.0;
                hopping_B_3 = 0.0;
            end
            ham = TopologicalInsulator_SSH.SSH_hamiltonian(hoppingA,hoppingB,sites,open,hopping_A_3,hopping_B_3);
            obj@TopologicalInsulator(ham,2);
            obj.hoppingA = hoppingA;
            obj.hoppingB = hoppingB;
            obj.hopping_A_3 = hopping_A_3;
            obj.hopping_B_3 = hopping_B_3;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
            ham_k = [[0, obj.hoppingA + (obj.hoppingB')*exp(-1i*k) + ...
                obj.hopping_A_3*exp(1i*k) + (obj.hopping_B_3')*exp(-2i*k)];[0,0]];
            ham_k = ham_k + ham_k';
        end
        
        function sps = BL_ground_spinors(obj,k_vals)
            sps = cell(1,numel(k_vals));
            for k_index = 1:numel(k_vals)
                k = k_vals(k_index);
                ham_k = obj.BL_k_hamiltonian(k);
                    assert(all(size(ham_k) == [2,2]));
                    assert(abs(ham_k(1,1)) < 1.e-10 && abs(ham_k(2,2)) < 1.e-10);
                sps{1,k_index} = sqrt(0.5)*[1;ham_k(2,1)/abs(ham_k(2,1))];
            end
        end
        
        
    end
    
    methods (Static)
        function ham = SSH_hamiltonian(hopping1, hopping2, sites,open,hopping_A_3,hopping_B_3)
            assert(mod(sites,2) == 0, 'Must have even number of sites');
            cells = sites/2;
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,[hopping1,hopping2],cells,open);
            ham3 = TopologicalInsulator.off_diagonal_matrix(3,[hopping_A_3,hopping_B_3],cells,open);
            ham = ham1 + ham3;
            ham = ham + ham';
        end
        
        function fig_handle = plot_current_entanglement(times,top_invars_1,top_invars_2,charges_1,charges_2,ent_1,ent_2,times_H_1,ent_H_1,times_H_2,ent_H_2)
            pt_to_inch = 72;
            width = 246; height = 140;
            fig_handle = figure('Name', 'Entanglement spectrum',...
                'Units','inches','Position',[4,3,width/pt_to_inch,height/pt_to_inch]);
            lw = 0.75;
            colblue = [0,0.4470,0.7410];
            colred = [0.8,0.2,0.05];
            colgrn = [0.1,0.7,0.2];
            colprp = [0.3,0.05,0.4];
            
            %SEARCH FOR DISCONTINUITIES
            d_1 = diff(top_invars_1);
            disc_1 = find(abs(d_1) > 0.5);
            if numel(disc_1) > 1
                warning('Ignored discontinuity');
                disc_1(2:end) = [];
            elseif numel(disc_1) == 0
                disc_1 = numel(top_invars_1);
            end
            
            d_1_m = diff(mod(charges_1,1));
            disc_1_m = find(abs(d_1_m) > 0.5);
            if numel(disc_1_m) > 1
                warning('Ignored discontinuity');
                disc_1_m(2:end) = [];
            elseif numel(disc_1_m) == 0
                disc_1_m = numel(charges_1);
            end
            
            
            
            hold on;
            h(1) = plot(times(1:disc_1),mod(top_invars_1(1:disc_1),1),'LineWidth',lw,'Color',colblue,'DisplayName','AIII');
            plot(times(1:disc_1_m),mod(charges_1(1:disc_1_m),1),'--','LineWidth',lw,'Color',colblue);
            plot(times((disc_1+1):end),mod(top_invars_1((disc_1+1):end),1),'LineWidth',lw,'Color',colblue);
            plot(times((disc_1_m+1):end),mod(charges_1((disc_1_m+1):end),1),'--','LineWidth',lw,'Color',colblue);
            h(2) = plot(times,mod(top_invars_2,1),'LineWidth',lw,'Color',colred,'DisplayName','BDI');
            plot(times,mod(charges_2,1),'--','LineWidth',lw,'Color',colred);
            %ylim([-0.8,0]);
            ylabel('CS$_1(t)$ mod 1; $Q(t)$ mod 1','interpreter','latex');
            hold off;
            xlabel('Time $t$ [units of $1/J_1$]','interpreter','latex');
            ax1 = gca;
            
            set(ax1,'TickLabelInterpreter','latex');
            set(ax1,'FontSize',8);
            set(ax1,'FontName','Times');
            
            set(ax1,'Units','normalized');
            pos1 = get(ax1,'Position');
            pos1(3) = 0.46;
            set(ax1,'Position',pos1);
            
            box on;
            l = legend(h,'Location','SouthEast');
            
            set(l,'Box','off');
            
            midmargin = 0.11;
            rightmargin = 0.01;
            
            ax2 = axes('Position',[pos1(1) + pos1(3) + midmargin,pos1(2),...
                1 - rightmargin - pos1(1) - pos1(3) - midmargin,pos1(4)/2 + 0.001]);
            ax3 = axes('Position',[pos1(1) + pos1(3) + midmargin,pos1(2) + pos1(4)/2,...
                1 - rightmargin - pos1(1) - pos1(3) - midmargin,pos1(4)/2]);
            
            
            
            hold on;
            plot(ax3,times,ent_1,'Color',colblue);
            plot(ax3,times,ent_2,'Color',colred);
            hold off;
            
            axes(ax2);
            
            hold on;
            plot(ax2,times_H_1,1.e3*(ent_H_1(2,:) - ent_H_1(1,:)),'Color',colgrn);
            plot(ax2,times_H_2,ent_H_2(2,:) - ent_H_2(1,:),'Color',colprp);
            hold off;
            
            
            
            etimemax = 3;
            xtic = [0,2];
            ylim2 = [0,0.15];
            ytic = [0,0.1];
            
            set(ax2,'XLim',[0,etimemax]);
            set(ax3,'XLim',[0,etimemax]);
            set(ax2,'XTick',xtic);
            set(ax3,'XTick',xtic);
            
            set(ax3,'Xticklabel',[]);
            
            set(ax2,'YLim',ylim2);
            set(ax2,'YTick',ytic);
            
            set(ax2,'TickLabelInterpreter','latex');
            set(ax2,'FontSize',8);
            set(ax2,'FontName','Times');
            set(ax3,'TickLabelInterpreter','latex');
            set(ax3,'FontSize',8);
            set(ax3,'FontName','Times');
            
            axes(ax2);
            
            xlabel('Time $t$ [units of $1/J_1$]','interpreter','latex');
            ylabel('Entanglement gap $\Delta_E$','interpreter','latex');
            
            box on;
            axes(ax3); box on;
            
            annotation('textbox','String','(a)','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',[0.2,0.2,0.2,0.1]);
            annotation('textbox','String','(b)','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',[0.2,0.2,0.2,0.1]);
            annotation('textbox','String','(c)','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',[0.2,0.2,0.2,0.1]);
            annotation('textbox','String','$\times 10^3$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','FontSize',8,'Position',[0.2,0.2,0.2,0.1]);
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

            phs = sum(sum(abs(corrmat + mat_phs)));
            trs = sum(sum(abs(corrmat - mat_trs)));
            chi = sum(sum(abs(corrmat + mat_chi)));
            
        end
    end
end

