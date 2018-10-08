classdef TopologicalInsulator_KM < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hopping
        so
        y_cells
    end
    
    methods
        function obj = TopologicalInsulator_KM(hopping,so,rashba,sublattice,kx,sites,open)
            ham = TopologicalInsulator_KM.KM_hamiltonian(hopping,so,rashba,sublattice,kx,sites,open);
            obj@TopologicalInsulator(ham,8);
            obj.hopping = hopping; obj.so = so;
            obj.y_cells = sites;
        end
        
        function ham_k = BL_k_hamiltonian(~,~)
            error('Not supported');
        end
        
        
    end
    
    methods(Static)
        %function KM_Hamiltonian: constructs a Kane-Mele Hamiltonian which
        %is fourier transformed in one direction only (along the point of
        %the hexagons). Spin labels are carried in the indices. The
        %ordering of indices is first cell number (y), then site index
        %(a-d), then spin (^ or v). i.e.
        %[a^, av, b^, bv, c^, cv, d^, dv, a(y+1)^, ...]
        %Note: sites here refers to the number of cells in the y direction
        function ham = KM_hamiltonian(hopping,so,r,sublattice,k,sites,open)
            %Notation: ^ = up spin, v = down spin. Letters a-d label the
            %4 sites within a cell.
            ek = exp(1i*k);
            %hoppings2 = [a^-b^, av-bv, b^-c^, bv-cv, c^-d^, cv-dv, ...
            %..., d^-a^(+1), dv-av (+1)];
            hoppings_2 = hopping*[1,1,1,1,1,1,exp(1i*k),exp(1i*k)];
            
            %hoppings4 = [a^-c^, av-cv, b^-d^, bv-dv, c^-a^(+1), cv-av(+1), ...
            %..., d^-b^(+1), dv-bv (+1)];
            hoppings_4 = 1i*so*[+1+ek',-1-ek',-1-ek',+1+ek',+1+ek,-1-ek,-1-ek,+1+ek];
            
            %hoppings6 = [a^-d^, av-dv, b^-a^(+1), bv-av(+1), c^-b^(+1), cv-bv(+1), ...
            %..., d^-c^(+1), dv-cv (+1)];
            hoppings_6 = hopping*[ek',ek',0,0,1,1,0,0];
            
            %hoppings8 = [a^-a^, av-av, b^-b^, bv-bv, c^-c^, cv-cv, ...
            %..., d^-d^, dv-dv] all (+1);
            hoppings_8 = 1i*so*[-1,+1,+1,-1,-1,+1,+1,-1];
            
            r3 = sqrt(3);
            
            %rahsba1 = [a^-av, av-b^, b^-bv, bv-c^, c^-cv, cv-d^, ...
            %..., d^-dv, dv-a^ (+1)];
            rashba_1 = 1i*r*[0,1i,0,r3-1i,0,1i,0,(r3-1i)*ek];
            
            %rahsba3 = [a^-bv, av-c^, b^-cv, bv-d^, c^-dv, cv-a^(+1), ...
            %..., d^-av(+1), dv-b^ (+1)];
            rashba_3 = 1i*r*[-1i,0,r3+1i,0,-1i,0,(r3+1i)*ek,0];
            
            %rahsba5 = [a^-cv, av-d^, b^-dv, bv-a^, c^-av, cv-b^(+1), ...
            %..., d^-bv(+1), dv-c^ (+1)];
            rashba_5 = 1i*r*[0,(r3+1i)*ek',0,0,0,r3+1i,0,0];
            
            %rahsba7 = [a^-dv, av-a^(+1), b^-av(+1), bv-b^(+1), c^-bv(+1), cv-c^(+1), ...
            %..., d^-cv(+1), dv-d^ (+1)];
            rashba_7 = 1i*r*[(r3-1i)*ek',0,0,0,r3-1i,0,0,0];
            
            on_site_pots = sublattice*[1,1,-1,-1,1,1,-1,-1];
            
            onsite = TopologicalInsulator.off_diagonal_matrix(0,on_site_pots,sites,open);
            hop2 = TopologicalInsulator.off_diagonal_matrix(2,hoppings_2,sites,open);
            hop4 = TopologicalInsulator.off_diagonal_matrix(4,hoppings_4,sites,open);
            hop6 = TopologicalInsulator.off_diagonal_matrix(6,hoppings_6,sites,open);
            hop8 = TopologicalInsulator.off_diagonal_matrix(8,hoppings_8,sites,open);
            
            hop1 = TopologicalInsulator.off_diagonal_matrix(1,rashba_1,sites,open);
            hop3 = TopologicalInsulator.off_diagonal_matrix(3,rashba_3,sites,open);
            hop5 = TopologicalInsulator.off_diagonal_matrix(5,rashba_5,sites,open);
            hop7 = TopologicalInsulator.off_diagonal_matrix(7,rashba_7,sites,open);
            
            
            ham = hop1 + hop2 + hop3 + hop4 + hop5 + hop6 + hop7 + hop8 + onsite/2;
            ham = ham + ham';
        end
        
        function trs = time_reversal_operator(y_cells)
            spin_exch = [[0, 1];[-1,0]];
            trs = kron(eye(y_cells*4),spin_exch);
        end
        
        function t = matrix_has_trs(mat_k,mat_km)
            assert(all(size(mat_k) == size(mat_km)),'Matrices must have same size');
            trs = TopologicalInsulator_KM.time_reversal_operator(size(mat_k,1)/8);
            t_conj = -trs * conj(mat_k) * trs;
            if sum(sum(abs(t_conj - mat_km))) < 1.e-6
                t = true;
            else
                disp(t_conj - mat_km);
                t = false;
            end
        end
        
        function fig_handle = entanglement_plot(kvals,spec_init_hal,spec_final_hal,spec_init_km,spec_final_km)
            width = 246; height = 200;
            fig_handle = figure('Name','Haldane and Kane Mele Entanglement Spectrum','Units','points',...
                'position',[300,200,width,height]);
            
            pos_1_km = spec_init_km(:,1) > 0; fp1 = find(pos_1_km); fp2 = find(~pos_1_km);
            pos_2_km = spec_final_km(:,1) > 0; fp3 = find(pos_2_km); fp4 = find(~pos_2_km);
            [~,ind_1_km] = min(spec_init_km(pos_1_km,1));
            [~,ind_2_km] = max(spec_init_km(~pos_1_km,1));
            [~,ind_3_km] = min(spec_final_km(pos_2_km,1));
            [~,ind_4_km] = max(spec_final_km(~pos_2_km,1));
            
            i1_km = fp1(ind_1_km);
            i2_km = fp2(ind_2_km);
            i3_km = fp3(ind_3_km);
            i4_km = fp4(ind_4_km);
            
            zind = find(abs(kvals) < 1.e-10);
            
            spc1_km = spec_init_km(i1_km,:);
            spc2_km = spec_init_km(i2_km,:);
            
            spc1_km(zind) = 0; spc2_km(zind) = 0;
            tmp = spc1_km(zind:end);
            spc1_km(zind:end) = spc2_km(zind:end);
            spc2_km(zind:end) = tmp;
            
             pos_1_hal = spec_init_hal(:,1) > 0; fp1 = find(pos_1_hal); fp2 = find(~pos_1_hal);
            pos_2_hal = spec_final_hal(:,1) > 0; fp3 = find(pos_2_hal); fp4 = find(~pos_2_hal);
            [~,ind_1_hal] = min(spec_init_hal(pos_1_hal,1));
            [~,ind_2_hal] = max(spec_init_hal(~pos_1_hal,1));
            [~,ind_3_hal] = min(spec_final_hal(pos_2_hal,1));
            [~,ind_4_hal] = max(spec_final_hal(~pos_2_hal,1));
            
            i1_hal = fp1(ind_1_hal);
            i2_hal = fp2(ind_2_hal);
            i3_hal = fp3(ind_3_hal);
            i4_hal = fp4(ind_4_hal);
            
            zind = find(abs(kvals) < 1.e-10);
            
            spc1_hal = spec_init_hal(i1_hal,:);
            spc2_hal = spec_init_hal(i2_hal,:);
            
            spc1_hal(zind) = 0; spc2_hal(zind) = 0;
            tmp = spc1_hal(zind:end);
            spc1_hal(zind:end) = spc2_hal(zind:end);
            spc2_hal(zind:end) = tmp;
            
            ent_func = @(x) (exp(x) + 1).^(-1);
            
            leftmargin = 15;
            bottommargin = 15;
            gap = 12;
            max_ent_km = 0.26;
            max_ent_hal = 0.52;
            max_k = max(kvals);
            x_limit = max_k*1.1;
            
            x_tick = 0.5;
            y_tick_km = 0.1;
            y_tick_hal = 0.2;
            y_limit_km = y_tick_km*1.2;
            y_limit_hal = y_tick_hal*1.2;
            b_label_width = 0.4;
            l_arrow_width = 0.3;
            
            col1 = [     0    0.4470    0.7410];
            col2 = [0.8500    0.3250    0.0980];
            
            fontsize = 8;
            
            plotwidth = (width - gap - leftmargin)/2;
            plotheight = (height - bottommargin - gap)/2;
            
            
            ax(1) = axes('units','points','position',[leftmargin,bottommargin + plotheight + gap,plotwidth,plotheight]);
            ax(2) = axes('units','points','position',[(leftmargin + gap + width)/2,bottommargin + plotheight + gap,plotwidth,plotheight]);
            
            ax(3) = axes('units','points','position',[leftmargin,bottommargin,plotwidth,plotheight]);
            ax(4) = axes('units','points','position',[(leftmargin + gap + width)/2,bottommargin,plotwidth,plotheight]);
            
            
            
            
            plot(ax(1),kvals,ent_func(spec_init_hal([i1_hal,i2_hal],:)),'color',col1);
            plot(ax(2),kvals,ent_func(spec_final_hal([i3_hal,i4_hal],:)),'color',col2);
            
            plot(ax(3),kvals,ent_func([spc1_km;spc2_km]),'color',col1);
            plot(ax(4),kvals,ent_func(spec_final_km([i3_km,i4_km],:)),'color',col2);
            
            for j = 1:4
                hold(ax(j),'on');
                plot(ax(j),[0,x_limit],[0.5,0.5],'--',...
                    'Color','black','LineWidth',0.25);
                plot(ax(j),[-x_limit*0.29,-x_limit],[0.5,0.5],'--',...
                    'Color','black','LineWidth',0.25);
                hold(ax(j),'off');
            end
            
            
            
            set(ax(1:2),'Ylim',[0.5-y_limit_hal,0.5+y_limit_hal]);
            set(ax(3:4),'Ylim',[0.5-y_limit_km,0.5+y_limit_km]);
            %set(ax,'XAxisLocation','origin');
            set(ax,'YAxisLocation','origin');
            set(ax,'Box','off');
            set(ax,'Xlim',[-x_limit,x_limit]);
            set(ax,'Xtick',[-x_tick,x_tick]);
            %set(ax1,'Xticklabel',[-max_k,max_k]);
            set(ax(1:2),'Ytick',[0.5-y_tick_hal,0.5,0.5+y_tick_hal]);
            set(ax(3:4),'Ytick',[0.5-y_tick_km,0.5,0.5+y_tick_km]);
            set(ax,'TickLabelInterpreter','latex');
            set(ax,'fontsize',fontsize);
            %xtickformat('%.2f');
    
            
            b_arrow_center = (width+leftmargin)/(2*width);
            l_arrow_center = (height + bottommargin)/(2*height);
            t_label_width = 0.2;
            l_label_space = (leftmargin*0.42/width);
            
            axis_shift_left = 0.064;
            axis_shift_up = 0.025;
            axis_shift_right = -0.01;
            axis_shift_down = 0.02;
            
            annotation('textbox','String','$k_x$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + plotwidth)/width) - axis_shift_left,(bottommargin)/height + axis_shift_up,...
                0.1,0.04],...
                'HorizontalAlignment','center','fontsize',fontsize);
            annotation('textbox','String','$k_x$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + 2*plotwidth + gap)/width) - axis_shift_left,(bottommargin)/height + axis_shift_up,...
                0.1,0.04],...
                'HorizontalAlignment','center','fontsize',fontsize);
            annotation('textbox','String','$k_x$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + plotwidth)/width) - axis_shift_left,(bottommargin + plotheight + gap)/height + axis_shift_up,...
                0.1,0.04],...
                'HorizontalAlignment','center','fontsize',fontsize);
            annotation('textbox','String','$k_x$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + 2*plotwidth + gap)/width) - axis_shift_left,(bottommargin + plotheight + gap)/height + axis_shift_up,...
                0.1,0.04],...
                'HorizontalAlignment','center','fontsize',fontsize);
            
            annotation('textbox','String','$\zeta_n$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + plotwidth/2)/width) + axis_shift_right,(bottommargin + plotheight)/height - axis_shift_down,...
                0.1,0.04],...
                'HorizontalAlignment','left','fontsize',fontsize);
            annotation('textbox','String','$\zeta_n$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + 3*plotwidth/2  + gap)/width) + axis_shift_right,(bottommargin + plotheight)/height - axis_shift_down,...
                0.1,0.04],...
                'HorizontalAlignment','left','fontsize',fontsize);
            annotation('textbox','String','$\zeta_n$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + plotwidth/2)/width) + axis_shift_right,(bottommargin + 2*plotheight + gap)/height - axis_shift_down,...
                0.1,0.04],...
                'HorizontalAlignment','left','fontsize',fontsize);
            annotation('textbox','String','$\zeta_n$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',...
                [((leftmargin + 3*plotwidth/2 + gap)/width) + axis_shift_right,(bottommargin + 2*plotheight + gap)/height - axis_shift_down,...
                0.1,0.04],...
                'HorizontalAlignment','left','fontsize',fontsize);
            
            arrow_len = 3;
            arrow_wid = 3;
            
            xshifts = [-0.002,-0.002,-0.001,-0.001];
            yshifts = [0.0005,0.000,0.0005,0.000];
            for j = 1:4
                axp = ax(j).Position;
                
                xe = (axp(1) + axp(3))/width;
                xs = xe - arrow_len/width;
                ye = (axp(2))/height + xshifts(j);
                
                annotation('arrow',[xs xe],[ye ye],'units','points',...
                    'headlength',arrow_len,'headwidth',arrow_wid);
                
                xe = (axp(1) + axp(3)/2)/width + yshifts(j);
                ye = (axp(2) + axp(4))/height;
                ys = ye - arrow_len/height;
                
                annotation('arrow',[xe xe],[ys ye],'units','points',...
                    'headlength',arrow_len,'headwidth',arrow_wid);
            end
            
            
%             annotation('textbox','String','Momentum $k_x/2\pi$','FitBoxToText','on',...
%                 'Linestyle','none','interpreter','latex','Position',[b_arrow_center - (b_label_width/2),0.08*bottommargin/height,b_label_width,0.65*bottommargin/height],...
%                 'HorizontalAlignment','center','fontsize',fontsize);
%             lab = ylabel('Entanglement energy $\epsilon_i$','interpreter','latex');
            t_hal = text(ax(1),-0.04,0.5,'Haldane','HorizontalAlignment','center','Rotation',90,'Units','normalized','Interpreter','latex');
            t_km = text(ax(3),-0.04,0.5,'Kane-Mele','HorizontalAlignment','center','Rotation',90,'Units','normalized','Interpreter','latex');
%             lab.Rotation = 90;
%             annotation('doublearrow',[b_arrow_center - (b_label_width/2),b_arrow_center + (b_label_width/2)],[0.02,0.02],...
%                 'Head1Length',4,'Head2Length',4,'Head1Width',3,'Head2Width',3);
%             annotation('doublearrow',[l_label_space,l_label_space],[l_arrow_center - (l_arrow_width/2),l_arrow_center + (l_arrow_width/2)],...
%                 'Head1Length',4,'Head2Length',4,'Head1Width',3,'Head2Width',3);
            
            annotation('textbox','String','$t=0$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',[(leftmargin + plotwidth/2)/width - t_label_width/2,0.5*bottommargin/height,t_label_width,0.5*bottommargin/height],...
                'HorizontalAlignment','center','fontsize',fontsize);
            annotation('textbox','String','$t=2J$','FitBoxToText','on',...
                'Linestyle','none','interpreter','latex','Position',[1 - (plotwidth/(2*width)) - t_label_width/2,0.5*bottommargin/height,t_label_width,0.5*bottommargin/height],...
                'HorizontalAlignment','center','fontsize',fontsize);
        end
        
    end
end

