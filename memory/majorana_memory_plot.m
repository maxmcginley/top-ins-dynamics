function fig_handle = majorana_memory_plot(fidelities,times)

width = 246;
height = 150;
leftmargin = 34;
rightmargin = 10;
bottommargin = 25;
topmargin = 6;

pos = [300,200,width,height];

fig_handle = figure('Name','Fidelity comparison','Units','points','Position',pos);

ax1 = axes('Units','points','Position',[leftmargin,bottommargin,width - leftmargin - rightmargin,...
    height - topmargin - bottommargin]);

leftcol = [0.9,0.97,1];
rightcol = [1,1,1];

% rectangle(ax1,'Position',[0,0.5,TIME_MAX,0.5],'FaceColor',leftcol,'EdgeColor',leftcol);
% rectangle(ax1,'Position',[TIME_MAX,0.5,100,0.5],'FaceColor',rightcol,'EdgeColor',rightcol);
% line(ax1,[TIME_MAX,TIME_MAX],[0.5,1],'LineStyle','--','Color','black');

set(ax1,'FontSize',8);
set(ax1,'FontName','Times');

hold(ax1,'on');
h(1) = plot(ax1,times,fidelities{1}/2,'DisplayName','(17a) -- DIII');
h(2) = plot(ax1,times,fidelities{2}/2,'DisplayName','(17b) -- BDI $\nu = 1$');
h(3) = plot(ax1,times,fidelities{3}/2,'DisplayName','(17c) -- BDI $\nu = 2$');
l = legend(h,'Location','SouthEast');
l.Interpreter = 'latex';

set(ax1,'TickLabelInterpreter','latex');
xlabel(ax1,'Time $t$ $[$Units of $J^{-1}]$','interpreter','latex');
ylabel(ax1,'Fidelity $\| \Gamma(\rho^+) - \Gamma(\rho^-)\|/2$','interpreter','latex');



%set(ax1,'YLim',[floor(min(output.fidelities{1}/2)*10)/10,1]);
set(ax1,'YLim',[-inf,1]);
set(ax1,'XLim',[0,max(times)]);

set(ax1,'Layer','top')

box(ax1,'on');

end