function maj_params = setup_insulators(maj_params)

if maj_params.num_insulators == 0
    maj_params.num_insulators = uint32(3);
end

ins_init = cell(1,maj_params.num_insulators);
hamiltonians = cell(1,maj_params.num_insulators);
spectra = cell(1,maj_params.num_insulators);

t = maj_params.system_params.t;
mu_init = maj_params.system_params.mu_init;
del_1 = maj_params.system_params.del_1;
del_2 = maj_params.system_params.del_2;
alpha = maj_params.system_params.alpha;
sites = maj_params.system_params.sites;
mu_patch_size = maj_params.system_params.mu_patch_size;

    ins_init{1} = TopologicalInsulator_DIII(t,mu_init,del_1,del_2,alpha,sites*4,true);
if maj_params.num_insulators >= 2
    ins_init{2} = TopologicalInsulator_DoubleKitaev(t,mu_init,del_1,sites*4,true);
end
if maj_params.num_insulators >= 3
    ins_init{3} = TopologicalInsulator_ChiralMaj(alpha,t,del_2,del_1,mu_init,sites*2,true);
end

mu_location_TRS = [[ones(1,mu_patch_size),zeros(1,sites - mu_patch_size)];...
                    [zeros(1,sites - mu_patch_size),ones(1,mu_patch_size)]];
mu_location_Double = [[ones(1,mu_patch_size),zeros(1,sites*2 - mu_patch_size)];... %Two copies of the chains
                       [zeros(1,sites),ones(1,mu_patch_size),zeros(1,sites - mu_patch_size)];...
                       [zeros(1,sites - mu_patch_size),ones(1,mu_patch_size),zeros(1,sites)];...
                       [zeros(1,sites),zeros(1,sites - mu_patch_size),ones(1,mu_patch_size)]];
mu_location_Chiral = [[ones(1,mu_patch_size),zeros(1,sites - mu_patch_size)];...
                    [zeros(1,sites - mu_patch_size),ones(1,mu_patch_size)]];

if maj_params.double_noise
    n_chan = [2,4,2];
    disp('Running with double noise');
else
    n_chan = [1,1,1];
end

%Defines where the final Hamiltonians are
hamiltonians{1} = cell(1,n_chan(1));    
for j = 1:(n_chan(1))
    hamiltonians{1}{2*j-1} = TopologicalInsulator_DIII(0,mu_location_TRS(j,:),0,0,0,sites*4,true,false).hamiltonian;
    hamiltonians{1}{2*j} = TopologicalInsulator_DIII(0,0,0,mu_location_TRS(j,:),0,sites*4,true,false).hamiltonian;
end

if maj_params.num_insulators >= 2
    hamiltonians{2} = cell(1,n_chan(2)); %Independent noise
for j = 1:(n_chan(2))
    hamiltonians{2}{2*j-1} = TopologicalInsulator_DoubleKitaev(0,mu_location_Double(j,:),0,sites*4,true,false).hamiltonian;
    hamiltonians{2}{2*j} = TopologicalInsulator_DoubleKitaev(0,0,mu_location_Double(j,:),sites*4,true,false).hamiltonian;
end
end

if maj_params.num_insulators >= 3
    hamiltonians{3} = cell(1,n_chan(3));
for j = 1:(n_chan(3))
    hamiltonians{3}{2*j-1} = TopologicalInsulator_ChiralMaj(0,0,0,0,mu_location_Chiral(j,:),sites*2,true,false).hamiltonian;
    hamiltonians{3}{2*j} = TopologicalInsulator_ChiralMaj(0,0,mu_location_Chiral(j,:),0,0,sites*2,true,false).hamiltonian;
end
end
    
    if maj_params.use_cutoff
        cutoff = maj_params.cutoffs;
    else
        cutoff = [];
    end

for j = 1:maj_params.num_insulators
spectra{j} = TimeEvolution_Noise.generate_poissonians(repmat(maj_params.spec_freq_widths,1,n_chan(j)),...
                                                    repmat(maj_params.spec_amps,1,n_chan(j)),...
                                                    repmat(cutoff,1,n_chan(j)));
end
 
maj_params.insulator_init = ins_init;
maj_params.hamiltonians = hamiltonians;
maj_params.spectra = spectra;

end

