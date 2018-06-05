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

%******************INPUT DATA*******************
sites = 72;
open = true;
mu_1 = 0.5;
delp = 1;
dels = 0.5;
alpha = 0.4;
TIME_STEP = 1;
TIME_MAX = 10;
POST_STEPS = 190;
quenches_per_step = 20;
mu_min = 0.0;
mu_max = 0.1;
num_mus = 25;
cell_size = 4;
hopping_range = 1;
%*********************************************

num_quenches = double(idivide(uint32(TIME_MAX*quenches_per_step),TIME_STEP,'ceil'));
num_times = (num_quenches / quenches_per_step) + POST_STEPS;
times = double((0:(num_times-1))).*TIME_STEP;

mu_2_vals = rand(1,num_mus*num_quenches)*(mu_max - mu_min)/2 + (mu_max + mu_min)/2;

cells = sites/cell_size;

ins_TRS_1 = TopologicalInsulator_DIII(mu_1,delp,dels,alpha,sites,open);
%ins_TRS_2 = TopologicalInsulator_DIII(mu,delp,dels_2,alpha,sites,open);

ins_Double_1 = TopologicalInsulator_DoubleKitaev(mu_1,delp,sites,open);

ins_TRS_1.spectrum
ins_Double_1.spectrum

majorana_limit = 0.05;

rxp = 0.5*ones(2);
rxm = 0.5*[[1,-1];[-1,1]];
%rxm = [[1,0];[0,0]];

state_TRS_minus = ins_TRS_1.coherent_majorana_qubit(majorana_limit,rxp);
state_TRS_plus = ins_TRS_1.coherent_majorana_qubit(majorana_limit,rxm);
state_Double_minus = ins_Double_1.coherent_majorana_qubit(majorana_limit,rxp);
state_Double_plus = ins_Double_1.coherent_majorana_qubit(majorana_limit,rxm);

%% Test symmetries

[init_trs,init_phs,init_chi] = TopologicalInsulator_DIII.test_symmetries(eye(size(state_TRS_minus)) - 2*state_TRS_minus);
%[init_trs,init_phs,init_chi] = TopologicalInsulator_DIII.test_symmetries(ins_TRS_1.hamiltonian);


double_phs = TopologicalInsulator_DoubleKitaev.test_phs(ins_Double_1.hamiltonian);
double_phs = TopologicalInsulator_DoubleKitaev.test_phs(eye(size(state_Double_minus)) - 2*state_Double_minus);

crit = 1.e-6;

fprintf('DIII: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);
fprintf('D: PHS = %d \n',abs(double_phs) < crit);

%% Time evolve

final_state_TRS_minus = zeros([size(state_TRS_minus),numel(times)]);
final_state_TRS_plus = zeros([size(state_TRS_plus),numel(times)]);
final_state_Double_minus = zeros([size(state_Double_minus),numel(times)]);
final_state_Double_plus = zeros([size(state_Double_plus),numel(times)]);

assert(times(1) == 0,'First time must be zero');

prog_handle = waitbar(0,'Time evolving...');

for mu_index = 1:num_mus
    
    final_state_TRS_minus_real = state_TRS_minus;
    final_state_TRS_plus_real = state_TRS_plus;
    final_state_Double_minus_real = state_Double_minus;
    final_state_Double_plus_real = state_Double_plus;

    waitbar((mu_index-1)/num_mus,prog_handle);
    
    for q_index = 1:num_quenches
        dis_index = (mu_index - 1)*num_quenches + q_index;
        
        ins_TRS_2 = TopologicalInsulator_DIII(mu_2_vals(dis_index),delp,dels,alpha,sites,open);
        ins_Double_2 = TopologicalInsulator_DoubleKitaev(mu_2_vals(dis_index),delp,sites,open);

        if sum(abs(ins_TRS_2.spectrum) < majorana_limit) < 4
            error('No Majoranas in final Hamiltonian');
        end
        
        prev_index = max(q_index-1,1);
        
        time = TIME_STEP*(q_index - prev_index);
        
        final_state_TRS_minus_real =  ...
            ins_TRS_2.time_evolve_correlation_matrix(final_state_TRS_minus_real,time);
        final_state_TRS_plus_real =  ...
            ins_TRS_2.time_evolve_correlation_matrix(final_state_TRS_plus_real,time);
        final_state_Double_minus_real = ...
            ins_Double_2.time_evolve_correlation_matrix(final_state_Double_minus_real,time);
        final_state_Double_plus_real =  ...
            ins_Double_2.time_evolve_correlation_matrix(final_state_Double_plus_real,time);
        
        if mod(q_index - 1,quenches_per_step) == 0
            t_index = (q_index - 1)/quenches_per_step + 1;
            final_state_TRS_minus(:,:,t_index) = ...
                final_state_TRS_minus(:,:,t_index) + (final_state_TRS_minus_real/num_mus);
            final_state_TRS_plus(:,:,t_index) = ...
                final_state_TRS_plus(:,:,t_index) + (final_state_TRS_plus_real/num_mus);
            final_state_Double_minus(:,:,t_index) = ...
                final_state_Double_minus(:,:,t_index) + (final_state_Double_minus_real/num_mus);
            final_state_Double_plus(:,:,t_index) = ...
                final_state_Double_plus(:,:,t_index) + (final_state_Double_plus_real/num_mus);
        end
    end
    
    pre_post_index = t_index;
    
    for p_step = 1:POST_STEPS
        t_index = pre_post_index + p_step;
        final_state_TRS_minus_real =  ...
            ins_TRS_2.time_evolve_correlation_matrix(final_state_TRS_minus_real,times(t_index) - times(t_index-1));
        final_state_TRS_plus_real =  ...
            ins_TRS_2.time_evolve_correlation_matrix(final_state_TRS_plus_real,times(t_index) - times(t_index-1));
        final_state_Double_minus_real =  ...
            ins_Double_2.time_evolve_correlation_matrix(final_state_Double_minus_real,times(t_index) - times(t_index-1));
        final_state_Double_plus_real =  ...
            ins_Double_2.time_evolve_correlation_matrix(final_state_Double_plus_real,times(t_index) - times(t_index-1));
        
        final_state_TRS_minus(:,:,t_index) = ...
            final_state_TRS_minus(:,:,t_index) + (final_state_TRS_minus_real/num_mus);
        final_state_TRS_plus(:,:,t_index) = ...
            final_state_TRS_plus(:,:,t_index) + (final_state_TRS_plus_real/num_mus);
        final_state_Double_minus(:,:,t_index) = ...
            final_state_Double_minus(:,:,t_index) + (final_state_Double_minus_real/num_mus);
        final_state_Double_plus(:,:,t_index) = ...
            final_state_Double_plus(:,:,t_index) + (final_state_Double_plus_real/num_mus);
    end
end

if ishandle(prog_handle)
    close(prog_handle);
end



TRS_fidelities = zeros(1,numel(times));
Double_fidelities = zeros(1,numel(times));

U_TRS = TopologicalInsulator_DIII.dirac_to_majorana_matrix(size(state_TRS_minus,1));
U_Double = TopologicalInsulator_DoubleKitaev.dirac_to_majorana_matrix(size(state_Double_minus,1));
%U_TRS_inv = inv(U_TRS);
%U_Double_inv = inv(U_Double);

for t_index = 1:numel(times)
    fid_matrix_TRS = conj(U_TRS) * (final_state_TRS_minus(:,:,t_index) - final_state_TRS_plus(:,:,t_index)) * U_TRS.';
    fid_matrix_Double = conj(U_Double) * (final_state_Double_minus(:,:,t_index) - final_state_Double_plus(:,:,t_index)) * U_Double.';
    
    %*********NON-LOCALITY*************
    left_sites = 1:(size(fid_matrix_Double,1)/2);
    right_sites = ((size(fid_matrix_Double,1)/2) + 1):(size(fid_matrix_Double,1));
    fid_matrix_Double_nonlocal = fid_matrix_Double;
    fid_matrix_Double_nonlocal(left_sites,left_sites) = 0;
    fid_matrix_Double_nonlocal(right_sites,right_sites) = 0;
    %******************************
    
    TRS_fidelities(1,t_index) = max(abs(eig(fid_matrix_TRS)));
    Double_fidelities(1,t_index) = max(abs(eig(fid_matrix_Double_nonlocal)));
end




%% Plotting

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

rectangle(ax1,'Position',[0,0.5,TIME_MAX,0.5],'FaceColor',leftcol,'EdgeColor',leftcol);
rectangle(ax1,'Position',[TIME_MAX,0.5,100,0.5],'FaceColor',rightcol,'EdgeColor',rightcol);
line(ax1,[TIME_MAX,TIME_MAX],[0.5,1],'LineStyle','--','Color','black');

hold(ax1,'on');
h(1) = plot(ax1,times,TRS_fidelities/2,'DisplayName','DIII');
h(2) = plot(ax1,times,Double_fidelities/2,'DisplayName','D');
l = legend(h,'Location','SouthWest');
l.Interpreter = 'latex';

set(ax1,'TickLabelInterpreter','latex');
xlabel(ax1,'Time $t$','interpreter','latex');
ylabel(ax1,'Fidelity $\| \Gamma(\rho^+) - \Gamma(\rho^-)\|/2$','interpreter','latex');

set(ax1,'YLim',[floor(min(TRS_fidelities/2)*10)/10,1]);
set(ax1,'XLim',[0,times(end) + TIME_STEP]);

set(ax1,'Layer','top')

box(ax1,'on');