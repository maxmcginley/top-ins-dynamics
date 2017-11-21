classdef TopologicalInsulator
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sites
        hamiltonian
        spectrum
        orbitals
    end
    
    
    
    methods
        function obj = TopologicalInsulator(hamiltonian)
            if size(hamiltonian,1) ~= size(hamiltonian,2)
                error('Hamiltonian is not square');
            end
            if any(any(abs(hamiltonian - hamiltonian') > 1.e-13))
                warning(['Sum of antihermitian part is ',...
                     num2str(sum(sum(abs(hamiltonian - hamiltonian'))))]);
                error('Hamiltonian is not Hermitian');
            end
            
            obj.sites = size(hamiltonian,1);
            obj.hamiltonian = (hamiltonian + hamiltonian')/2;
            [Evec, Eval] = eig(obj.hamiltonian);
            obj.orbitals = Evec';
            obj.spectrum = real(diag(Eval));
        end
        
        function fig_handle = plot_orbitals_in_energy_range(obj,lowerbound,upperbound)
            orb_choices = find(obj.spectrum > lowerbound & obj.spectrum < upperbound);
            if sum(orb_choices) == 0
                fig_handle = [];
                sprintf('No orbitals within range [%f , %f]',lowerbound,upperbound);
                return;
            end
            fig_handle = figure('Name','Orbital wavefunction');
            plot(1:obj.sites, abs(obj.orbitals(orb_choices,:).'));
            xlabel('Site index');
            ylabel('Amplitude');
            xlim([1, obj.sites]);
        end
        
        
        
        function vec = zero_mode_wavefunction(obj,lowerbound,upperbound)
            orb_choices = find(obj.spectrum > lowerbound & obj.spectrum < upperbound, 1);
            if isempty(orb_choices)
                vec = NaN(1,obj.sites);
                sprintf('No orbitals within range [%f , %f]',lowerbound,upperbound);
                return;
            end
            vec = obj.orbitals(orb_choices,:);
        end
        
        %Defined at <Psi | c_j e^{+iHt} c_j^\dagger c_j e^{-iHt} c_j^\dagger | Psi>
%         where c_j is some linear combination of fermion operators
        %site and |Psi> is the slater determinant defined by the
        %correlation matrix corrmat
        function val = noneq_qubit_function(obj,corrmat,site,time)
            phase_matrix = diag(exp(1i*time*obj.spectrum));
            wavefunction_evol = obj.orbitals' * phase_matrix * obj.orbitals;
            
            corrmat_00 = eye(obj.sites) - corrmat; % Gives <c_j c_k^\dagger >
            corrmat_0t = (wavefunction_evol*(eye(obj.sites) - corrmat)).';
            corrmat_t0 = conj(wavefunction_evol)*(eye(obj.sites) - corrmat);
            corrmat_tt = wavefunction_evol * corrmat * wavefunction_evol';
            
            val = corrmat_0t(site,site)*corrmat_t0(site,site) ...
                + corrmat_tt(site,site)*corrmat_00(site,site);
        end
        
        %If chem_potential is NaN, fills exactly half of the states
        function mat = half_filled_correlation_matrix(obj,chem_potential)
            if isnan(chem_potential)
                if mod(obj.sites,2) == 0
                    mat = obj.correlation_matrix([ones(1,obj.sites/2),zeros(1,obj.sites/2)]);
                else
                    mat = obj.correlation_matrix([ones(1,(obj.sites-1)/2),zeros(1,(obj.sites+1)/2)]);
                end
            else
                mat = obj.correlation_matrix(obj.spectrum < chem_potential);
            end
        end
        
        function ham_fic = time_evolve_hamiltonian(obj,ham0,time)
            phases = exp(-1i*time*obj.spectrum);
            t_evol = obj.orbitals' * diag(phases) * obj.orbitals;
            ham_fic = t_evol * ham0 * t_evol';
        end
        
        function mat = correlation_matrix(obj,occs)
            mat = obj.orbitals' * diag(occs) * obj.orbitals;
        end
        
        function mat = time_evolve_correlation_matrix(obj,init_cor_mat,time)
            phase_matrix = diag(exp(1i*time*obj.spectrum));
            mat = obj.orbitals' * phase_matrix * obj.orbitals ...
                * init_cor_mat * obj.orbitals' * conj(phase_matrix) * obj.orbitals;
        end
        
        function vec = time_evolve_creation_operator(obj,init_wavefunction,time)
            phase_matrix = diag(exp(1i*time*obj.spectrum));
            vec = obj.orbitals' * phase_matrix * obj.orbitals ...
                * init_wavefunction;
        end
        
        
    end
    
    methods(Static)
        function ham = SSH_hamiltonian(hopping1, hopping2, sites,open)
            if mod(sites,2) == 0
                ham = diag([repmat([hopping1, hopping2],1,(sites-2)/2), hopping1], 1);
                if ~open
                    ham(sites,1) = hopping2;
                end
                
                ham = ham + ham';
            else
                if ~open
                    error('Must have even number of sites for periodic system');
                end
                ham = diag(repmat([hopping1, hopping2],1,(sites-1)/2), 1);
                ham = ham + ham';
            end
        end
        
        function mat = off_diagonal_matrix(off_index,vals,cells,open)
            if size(vals,1) ~= 1
                error('Values should be given as row matrix');
            end
            if off_index == 0
                mat = diag(repmat(vals,1,cells));
                return;
            end
            cell_size = size(vals,2);
            mat_size = cell_size*cells;
            corner_size = abs(off_index);
            principal_cells = idivide(mat_size - corner_size,uint32(cell_size),'floor');
            resid = mat_size - corner_size - principal_cells*cell_size;
            principal_vec = [repmat(vals,1,principal_cells),vals(1:(resid))];
            
            mat = diag(principal_vec,off_index);
            if ~open
                rem_vec = [vals((resid+1):cell_size),repmat(vals,1,cells - principal_cells - 1)];
                if off_index > 0
                    mat((mat_size+1-corner_size):mat_size,1:corner_size) = diag(rem_vec);
                else
                    mat(1:corner_size,(mat_size+1-corner_size):mat_size) = diag(rem_vec);
                end
            end
        end
        
        

        function ham = kitaev_chain_hamiltonian(mu, t, delta, sites, open)
            ham = diag(repmat(0.5*[mu, -mu],1,sites));
            ham = ham + diag(repmat(0.5*[t, -conj(t)],1,sites-1),2);
            ham = ham + diag([repmat(0.5*[0, delta],1,sites-1), 0],1);
            ham = ham + diag([repmat(0.5*[-conj(delta), 0],1,sites-2), -0.5*conj(delta)],3);
            
            if ~open
                ham(2*sites - 1,1) = 0.5*t;
                ham(2*sites,2) = -0.5*conj(t);
                ham(2*sites,1) = 0.5*delta;
                ham(2*sites-1,2) = -0.5*conj(delta);
            end
            
            ham = ham + ham';
        end
        
        function dens = local_chern_densities(corrmat,cell_size)
            sites = size(corrmat,1);
            if mod(sites,cell_size) ~= 0
                error('Correlation matrix size incommensurate with cell size');
            end
            
            cells = sites/cell_size;
            dens = zeros(1,cells);
            pos_values = exp(2i*pi*kron(0:(cells-1),ones(1,cell_size))/cells);
            %pos_values = exp(((0:(sites-1)))*2i*pi/sites);
            prodmat = logm((corrmat*diag(pos_values)*corrmat)+ eye(size(corrmat)) - corrmat);
            imag(diag(prodmat))*cells/(2*pi)
            sum(imag(diag(prodmat))/(2*pi))
            pd = diag(prodmat);
            for j = 1:(cells)
                cell_pos_values = (j-1)*cell_size + (1:(cell_size));
                  dens(j) = sum(imag(pd(cell_pos_values))*cells/(2*pi));
            end
            dens = mod(dens,1);
        end
%         
%         function bloch_curr = bloch_current_vals(ham,cell_size,occ_bands)
%             [evec,eval] = eig(ham);
%             sites = size(ham,1);
%             choices = 1:sites;
%             ens = diag(eval);
%             empties = ens > 0;
%             choices(empties) = [];
%             if numel(choices) ~= (sites/cell_size)*occ_bands
%                 error('Wrong number of filled states for occupied bands');
%             end
%             k_vals = mod(angle(evec(1+cell_size,choices)) - angle(evec(1,choices)),2*pi);
%             [sorted_ks, ind] = sort(k_vals);
%             bloch_curr = 0;
%             
%             bloch_states_m = ((sites/cell_size)-1)*cell_size + (1:occ_bands);
%             [~,ind_k_m] = sort(ens(bloch_states_m));
%             
%             for k_index = 1:(sites/cell_size)
%                 bloch_states = (k_index-1)*cell_size + (1:occ_bands);
%                 [~,ind_k] = sort(ens(bloch_states));
%                 for b_index = 1:occ_bands
%                     unk = evec(:,ind(bloch_states(ind_k(b_index))));
%                     proj_k = unk * unk';
%                     unkm = evec(:,ind(bloch_states_m(ind_k_m(b_index))));
%                     proj_k_m = unkm * unkm';
%                     deluk = (proj_k - proj_k_m)/(sorted_ks(bloch_states(1)) - sorted_ks(bloch_states_m(1)));
%                     del_ens = ens(choices(
%                     bloch_curr = bloch_curr - deluk' * ham * uk - uk' * ham * deluk';
%                 end
%                 ind_k_m = ind_k;
%                 bloch_states_m = bloch_states;
%             end
%         end
        
        function spec = entanglement_spectrum_from_correlation_matrix(corrmat, site1, site2)
            redcorrmat = corrmat(site1:site2,site1:site2);
            evals = eig(redcorrmat);
            spec = log(1 - evals) - log(evals);
        end
        
        function ent = entanglement_entropy_from_correlation_matrix(corrmat, site1, site2)
            redcorrmat = corrmat(site1:site2,site1:site2);
            evals = eig(redcorrmat);
            evals(abs(evals) < 1.e-9 | abs(evals - 1) < 1.e-9) = [];
            H2 = @(x) -x*log(x) - (1-x)*log(1-x);
            ent = sum(arrayfun(H2,evals));
        end
        
        function v = two_particle_expectation(corrmat,creation_vec,annihilation_vec)
            if size(creation_vec,1) ~= 1 || size(annihilation_vec,1) ~= 1
                error('Wavefunctions to be provided as row vectors');
            end
            v = creation_vec * corrmat * annihilation_vec.';
        end
        
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
            mat_phs = TopologicalInsulator.nambu_operator(1,sites) * conj(corrmat) * TopologicalInsulator.nambu_operator(1,sites);
            mat_trs = conj(corrmat);
            mat_chi = TopologicalInsulator.nambu_operator(3,sites) * corrmat * TopologicalInsulator.nambu_operator(3,sites);

            phs = sum(sum(abs((eye(sites) - 2*corrmat) + (eye(sites) - 2*mat_phs))));
            trs = sum(sum(abs((eye(sites) - 2*corrmat) - (eye(sites) - 2*mat_trs))));
            chi = sum(sum(abs((eye(sites) - 2*corrmat) + (eye(sites) - 2*mat_chi))));
            
        end
        
        function xs_sorted = wannier_centres(corrmat,prev_centres)
            sites = size(corrmat,1);
            pos_values = exp(((1:(sites)))*2i*pi/sites);
            proj_pos = corrmat*diag(pos_values)*corrmat;
            centres = eig(proj_pos);
            xs = angle(centres(abs(centres) > 1.e-6))*sites/(2*pi);
            
            if any(isnan(prev_centres))
                xs_sorted = sort(xs);
            else
                xs_sorted = zeros(numel(xs),1);
                xs_rem = xs;
                ind_rem  = 1:numel(xs);
                for j = 1:numel(xs)
                    mod_xs = [xs_rem, xs_rem + sites, xs_rem - sites, xs_rem + 2*sites, xs_rem - 2*sites];
                    [diffs,mods] = min(abs(mod_xs - prev_centres(j)),[],2);
                    %diffs = abs(exp(2i*pi*xs_rem) - exp(2i*pi*prev_centres(j)));
                    [~,ind] = min(diffs);
                    xs_sorted(j) = mod_xs(ind,mods(ind));
                    xs_rem(ind) = []; ind_rem(ind) = [];
                end
            end
        end
        
        % Current operator defined from the rate of change of a^dagger_site
        % a_site
%         function c_matrix = current_operator_from_hamiltonian(ham,site,open)
%             h_sites = size(ham,1);
%             sources = (site - uint32(floor(h_sites))):site;
%             if open
%                 sources(sources < 1) = [];
%             else
%                 sources = mod(sources - 1,h_sites) + 1;
%             end
%             c_matrix = zeros(h_sites);
%             for j = 1:numel(sources)
%                 source = sources(j);
%                 targets = (site+1):(source + uint32(floor(h_sites)));
%                 if open
%                     targets(targets > h_sites) = [];
%                 else
%                     targets = mod(targets - 1,h_sites) + 1;
%                 end
%                 c_matrix(targets,source) = -1i*ham(targets,source);
%             end
%             c_matrix = c_matrix + c_matrix';
%         end
        function c_matrix = current_operator_from_hamiltonian(ham,site,max_range,open)
            c_matrix = zeros(size(ham));
            sites = size(ham,1);
            if site <= max_range
                warning('Current site near periodic boundary');
            end
            for src = (site-max_range):(site-1)
                c_matrix(src,site:(src+max_range)) = -1i*ham(src,site:(src+max_range));
            end
            c_matrix = 1i*(c_matrix - c_matrix');
        end
        
        function c_matrix = current_operator_from_hamiltonian_test(ham,site,max_range,open)
            if site <= max_range
                warning('Current site within hopping range of boundary');
            end
            sites = size(ham,1);
            occs = [ones(1,site), zeros(1,sites-site)];
            comm = 1i*(diag(occs) * ham - ham*diag(occs));
            h_mask = zeros(size(comm));
            h_mask((site-max_range):(site+max_range),...
                (site-max_range):(site+max_range)) = ones(2*max_range + 1);
            c_matrix = comm .* h_mask;
        end
        
        
        
        function t = validate_current_operator(ham,site1,site2,max_range,open)
            c_mat_1 = TopologicalInsulator.current_operator_from_hamiltonian(ham,site1,max_range,open);
            c_mat_2 = TopologicalInsulator.current_operator_from_hamiltonian(ham,site2,max_range,open);
            occs = zeros(1,size(ham,1));
            occs(site1:(site2-1)) = ones(1,site2-site1);
            n_op = diag(occs);
            t = sum(sum(abs(c_mat_2 - c_mat_1 + 1i*n_op*ham - 1i*ham*n_op)));
        end
        
        
        function fig_handle = correlation_length_plot(corrmat,open)
            if size(corrmat,1) ~= size(corrmat,2)
                error('Correlation matrix is not square');
            end
            sites = size(corrmat,1);
            
            xdata = NaN(1,(sites*(sites+1))/2);
            ydata = NaN(1,(sites*(sites+1))/2);
            pointer = 1;
            for j = 1:sites
                for k = 1:j
                    corval = corrmat(j,k);
                    if ~open
                        xdata(pointer) = min([abs(j - k),abs(j - k - sites),abs(j-k+sites)]);
                    else
                        xdata(pointer) = j - k;
                    end
                    ydata(pointer) = abs(corval);
                    pointer = pointer + 1;
                end
            end
            
            fig_handle = figure('Name','Correlation length');
            semilogy(xdata,ydata,'o');
            xlabel('Distance');
            ylabel('Correlation');
        end
        
        
    end
end

