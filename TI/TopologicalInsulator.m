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
            t_evol = expm(-1i*time*obj.hamiltonian);
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
        
        % Current operator defined from the rate of change of a^dagger_site
        % a_site
        function c_matrix = current_operator_from_hamiltonian(ham,site,open)
            h_sites = size(ham,1);
            sources = (site - uint32(floor(h_sites))):site;
            if open
                sources(sources < 1) = [];
            else
                sources = mod(sources - 1,h_sites) + 1;
            end
            c_matrix = zeros(h_sites);
            for j = 1:numel(sources)
                source = sources(j);
                targets = (site+1):(source + uint32(floor(h_sites)));
                if open
                    targets(targets > h_sites) = [];
                else
                    targets = mod(targets - 1,h_sites) + 1;
                end
                c_matrix(targets,source) = -1i*ham(targets,source);
            end
            c_matrix = c_matrix + c_matrix';
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

