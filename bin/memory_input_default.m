maj_params = MajoranaMemory_Params();

maj_params.timestep = 0.008;
maj_params.num_steps = 150000;
maj_params.steps_per_measure = 1000;

maj_params.max_exp = 2;

maj_params.reals = 30;

maj_params.run_parallel = true;

system_params = struct();

system_params.sites = 16;
system_params.t = 1;
system_params.mu_init = 0.25;
system_params.del_1 = 1;
system_params.del_2 = 0.3;
system_params.alpha = 0.2;
system_params.mu_patch_size = 2;

system_params.double_noise = true;

maj_params.system_params = system_params;

maj_params.num_insulators = 1;

maj_params.spec_freq_widths = [0.01,0.01];

amplitudes = [[0.5,0.5];[0.6,0.6];[0.7,0.7];[0.8,0.8]];
%amplitudes = [[0.1,0.3];[0.6,0.6]];
maj_params.spec_amps = [0.1,0.1];

%maj_params.spec_amps = 

maj_params.use_cutoff = false;
maj_params.cutoffs = [0.2,0.2];
