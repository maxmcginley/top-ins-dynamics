function fidelities = compute_majorana_fidelities(final_state_minus,final_state_plus)

num_insulators = numel(final_state_minus);
num_vals = size(final_state_minus{1},3);

fidelities = cell(1,num_insulators);
for i = 1:num_insulators
    fidelities{i} = zeros(1,num_vals);
end

for t_index = 1:num_vals
%     fid_matrix_TRS = conj(U_TRS) * (final_state_TRS_minus(:,:,t_index) - final_state_TRS_plus(:,:,t_index)) * U_TRS.'/2;
%     fid_matrix_Double = conj(U_Double) * (final_state_Double_minus(:,:,t_index) - final_state_Double_plus(:,:,t_index)) * U_Double.'/2;
    
    

    for i = 1:num_insulators
        fid_matrix = final_state_minus{i}(:,:,t_index) - final_state_plus{i}(:,:,t_index);
        fidelities{i}(t_index) = sqrt(trace(fid_matrix * fid_matrix'));
    end
    
%     %*********NON-LOCALITY*************
%     left_sites_Double = [1:(size(fid_matrix_Double,1)/4), ((size(fid_matrix_Double,1)/2) + 1):(3*size(fid_matrix_Double,1)/4)];
%     right_sites_Double = [((size(fid_matrix_Double,1)/4)+1):((size(fid_matrix_Double,1)/2)), ((3*size(fid_matrix_Double,1)/4) + 1):(size(fid_matrix_Double,1))];
%     left_sites_TRS = 1:(size(fid_matrix_TRS,1)/2);
%     right_sites_TRS = ((size(fid_matrix_TRS,1)/2) + 1):(size(fid_matrix_TRS,1));
%     
%     fid_matrix_TRS_nonlocal = fid_matrix_TRS;
%     fid_matrix_TRS_nonlocal(left_sites_TRS,left_sites_TRS) = 0;
%     fid_matrix_TRS_nonlocal(right_sites_TRS,right_sites_TRS) = 0;
%     fid_matrix_Double_nonlocal = fid_matrix_Double;
%     fid_matrix_Double_nonlocal(left_sites_Double,left_sites_Double) = 0;
%     fid_matrix_Double_nonlocal(right_sites_Double,right_sites_Double) = 0;
%     %******************************
%     
%     TRS_fidelities(1,t_index) = max(abs(eig(fid_matrix_TRS_nonlocal)));
%     Double_fidelities(1,t_index) = max(abs(eig(fid_matrix_Double_nonlocal)));
end
end

