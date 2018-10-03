if exist('figure_handles','var') 
    for j = 1:numel(figure_handles)
        if ishandle(figure_handles{j})
            close(figure_handles{j});
        end
    end
    clear('figure_handles');
end

clc;
clear;

figure_handles = cell(1,1);

addpath(fullfile(pwd,'..','TI'));
addpath(fullfile(pwd,'..','TE'));
addpath(fullfile(pwd,'..','memory'));

filename = 'majorana_memory_out.mat';

%******************INPUT DATA*******************
maj_params = MajoranaMemory_Params();

maj_params.timestep = 0.025;
maj_params.num_steps = 200;
maj_params.steps_per_measure = 10;

maj_params.max_exp = 2;

maj_params.reals = 10;

maj_params.run_parallel = false;

system_params = struct();

system_params.sites = 16;
system_params.t = 1;
system_params.mu_init = 0.25;
system_params.del_1 = 1;
system_params.del_2 = 0.3;
system_params.alpha = 0.2;
system_params.mu_patch_size = 2;

maj_params.system_params = system_params;

maj_params.spec_freq_widths = [200,200];
maj_params.spec_amps = [0.25,0.15];

RUN = false;
SAVE = true;
%*********************************************

if RUN
    maj_params = setup_insulators(maj_params);    

    [final_state_minus,final_state_plus] = calculate_majorana_evolution (maj_params);

    fidelities = compute_majorana_fidelities(final_state_minus,final_state_plus);
end
%% Saving data

if RUN && SAVE
    output = struct('final_state_minus',{final_state_minus},'final_state_plus',{final_state_plus},...
        'maj_params',maj_params,'fidelities',{fidelities});
    save(filename,'-struct','output');
end

%% Plotting

if ~RUN
output = load(filename);


width = 246;
height = 150;
leftmargin = 34;
rightmargin = 10;
bottommargin = 25;
topmargin = 6;

pos = [300,200,width,height];

figure_handles{end+1} = figure('Name','Fidelity comparison','Units','points','Position',pos);

ax1 = axes('Units','points','Position',[leftmargin,bottommargin,width - leftmargin - rightmargin,...
    height - topmargin - bottommargin]);

leftcol = [0.9,0.97,1];
rightcol = [1,1,1];

% rectangle(ax1,'Position',[0,0.5,TIME_MAX,0.5],'FaceColor',leftcol,'EdgeColor',leftcol);
% rectangle(ax1,'Position',[TIME_MAX,0.5,100,0.5],'FaceColor',rightcol,'EdgeColor',rightcol);
% line(ax1,[TIME_MAX,TIME_MAX],[0.5,1],'LineStyle','--','Color','black');

hold(ax1,'on');
h(1) = plot(ax1,output.maj_params.data_times(),output.fidelities{1}/2,'DisplayName','DIII');
h(2) = plot(ax1,output.maj_params.data_times(),output.fidelities{2}/2,'DisplayName','D');
h(3) = plot(ax1,output.maj_params.data_times(),output.fidelities{3}/2,'DisplayName','BDI');
l = legend(h,'Location','SouthWest');
l.Interpreter = 'latex';

set(ax1,'TickLabelInterpreter','latex');
xlabel(ax1,'Time $t$','interpreter','latex');
ylabel(ax1,'Fidelity $\| \Gamma(\rho^+) - \Gamma(\rho^-)\|/2$','interpreter','latex');

%set(ax1,'YLim',[floor(min(output.fidelities{1}/2)*10)/10,1]);
set(ax1,'YLim',[-inf,1]);
set(ax1,'XLim',[0,max(output.maj_params.data_times())]);

set(ax1,'Layer','top')

box(ax1,'on');

end