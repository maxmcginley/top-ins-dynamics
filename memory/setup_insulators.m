function maj_params = setup_insulators(maj_params)

maj_params.num_insulators = uint32(3);
maj_params.num_channels = uint32(2);

ins_init = cell(1,maj_params.num_insulators);
hamiltonians = cell(maj_params.num_insulators,maj_params.num_channels);

t = maj_params.system_params.t;
mu_init = maj_params.system_params.mu_init;
del_1 = maj_params.system_params.del_1;
del_2 = maj_params.system_params.del_2;
alpha = maj_params.system_params.alpha;
sites = maj_params.system_params.sites;
mu_patch_size = maj_params.system_params.mu_patch_size;

ins_init{1} = TopologicalInsulator_DIII(t,mu_init,del_1,del_2,alpha,sites*4,true);
ins_init{2} = TopologicalInsulator_DoubleKitaev(t,mu_init,del_1,sites*4,true);
ins_init{3} = TopologicalInsulator_ChiralMaj(alpha,t,del_2,del_1,mu_init,sites*2,true);

mu_location_TRS = [ones(1,mu_patch_size),zeros(1,sites - mu_patch_size)];
mu_location_Double = [ones(1,mu_patch_size),zeros(1,sites*2 - mu_patch_size)]; %Two copies of the chains
mu_location_Chiral = [ones(1,mu_patch_size),zeros(1,sites - mu_patch_size)];

%Defines where the final Hamiltonians are
    hamiltonians{1,1} = TopologicalInsulator_DIII(0,mu_location_TRS,0,0,0,sites*4,true,false).hamiltonian;
    hamiltonians{2,1} = TopologicalInsulator_DoubleKitaev(0,mu_location_Double,0,sites*4,true,false).hamiltonian;
    hamiltonians{3,1} = TopologicalInsulator_ChiralMaj(0,0,0,0,mu_location_Chiral,sites*2,true,false).hamiltonian;
    
    hamiltonians{1,2} = TopologicalInsulator_DIII(0,0,0,mu_location_TRS,0,sites*4,true,false).hamiltonian;
    hamiltonians{2,2} = TopologicalInsulator_DoubleKitaev(0,0,mu_location_Double,sites*4,true,false).hamiltonian;
    hamiltonians{3,2} = TopologicalInsulator_ChiralMaj(0,0,mu_location_Chiral,0,0,sites*2,true,false).hamiltonian;

maj_params.insulator_init = ins_init;
maj_params.hamiltonians = hamiltonians;

end

