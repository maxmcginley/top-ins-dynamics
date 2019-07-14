function test_symmetries(states_in, num_ins)
%% Test symmetries

crit = 1.e-5;

[init_trs,init_phs,init_chi] = TopologicalInsulator_DIII.test_symmetries(eye(size(states_in{1})) - 2*states_in{1});
%[init_trs,init_phs,init_chi] = TopologicalInsulator_DIII.test_symmetries(ins_2{1}.hamiltonian);

fprintf('DIII: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs) < crit,abs(init_trs) < crit,abs(init_chi) < crit);


if num_ins >= 2
    % double_phs = TopologicalInsulator_DoubleKitaev.test_phs(ins_2{2}.hamiltonian);
    double_phs = TopologicalInsulator_DoubleKitaev.test_phs(eye(size(states_in{2})) - 2*states_in{2});
    fprintf('D: PHS = %d \n',abs(double_phs) < crit);
end

if num_ins >= 3
    [init_trs_BDI,init_phs_BDI,init_chi_BDI] = TopologicalInsulator_ChiralMaj.test_symmetries(eye(size(states_in{3})) - 2*states_in{3});
    % %[init_trs_BDI,init_phs_BDI,init_chi_BDI] = TopologicalInsulator_ChiralMaj.test_symmetries(ins_1{3}.hamiltonian);
    fprintf('BDI: PHS = %d ; TRS = %d ; CHI = %d \n',abs(init_phs_BDI) < crit,abs(init_trs_BDI) < crit,abs(init_chi_BDI) < crit);
end

end