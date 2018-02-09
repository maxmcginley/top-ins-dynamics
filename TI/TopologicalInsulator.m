classdef (Abstract) TopologicalInsulator
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sites
        hamiltonian
        spectrum
        orbitals
        cell_size
    end
    
    
    
    methods
        %*************CONSTRUCTOR**********************
        function obj = TopologicalInsulator(hamiltonian, cell_size)
            if size(hamiltonian,1) ~= size(hamiltonian,2)
                error('Hamiltonian is not square');
            end
            assert(mod(size(hamiltonian,1),cell_size) == 0,...
                'Sites must be a multiple of cell size');
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
            obj.cell_size = cell_size;
        end
        
        function fig_handle = plot_orbitals_in_energy_range(obj,lowerbound,upperbound)
            orb_choices = find(obj.spectrum > lowerbound & obj.spectrum < upperbound);
            if sum(orb_choices) == 0
                fig_handle = [];
                sprintf('No orbitals within range [%f , %f]',lowerbound,upperbound);
                return;
            end
            fig_handle = obj.plot_orbitals_by_index(orb_choices);
        end
        
        function fig_handle = plot_orbitals_by_index(obj,orb_choices)
            fig_handle = figure('Name','Orbital wavefunction');
            plot(1:obj.sites, abs(obj.orbitals(orb_choices,:).'));
            xlabel('Site index');
            ylabel('Amplitude');
            xlim([1, obj.sites]);
        end
        
        %**************CORRELATION MATRIX***********************
        
        function mat = correlation_matrix(obj,occs)
            mat = obj.orbitals' * diag(occs) * obj.orbitals;
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
        
        
        
        %***************TIME EVOLUTION**************************
        
        function ham_fic = time_evolve_hamiltonian(obj,ham0,time)
            phases = exp(-1i*time*obj.spectrum);
            t_evol = obj.orbitals' * diag(phases) * obj.orbitals;
            ham_fic = t_evol * ham0 * t_evol';
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
        
        function nus = BL_topological_invariant(obj,init_spinors,times,k_vals)
            assert(iscell(init_spinors),'Initial spinors must be provided as a cell');
            assert(numel(init_spinors) == numel(k_vals),...
                'Number of k values must match number of spinors');
            occ_bands = size(init_spinors{1},2);
            assert(size(init_spinors{1},1) == obj.cell_size,'Spinor size must equal cell size');
            
            bloch_vectors = zeros(obj.cell_size,occ_bands,numel(k_vals),numel(times));
            for j = 1:numel(k_vals)
                k = k_vals(j);
                hamilt_k = obj.BL_k_hamiltonian(k);
                for t_index = 1:numel(times)
                    time = times(t_index);
                    evol = expm(-1i*hamilt_k*time);
                    bloch_vectors(:,:,j,t_index) = evol * init_spinors{j};
                end
            end
            
            nus = zeros(1,numel(times));
            
            for t_index = 1:numel(times)
                nus(1,t_index) = TopologicalInsulator.BL_wilson_loops(...
                    bloch_vectors(:,:,:,t_index));
            end
        end
    end
    
    methods (Abstract)
        ham_k = BL_k_hamiltonian(obj,k)
    end
    
    methods(Static)
        
        function sps = BL_constant_spinor(spinor,k_vals)
            sps = cell(1,numel(k_vals));
            sps(:) = {spinor};
        end
        
        function wilson_loop = BL_wilson_loops(spins)
            cell_size = size(spins,1);
            occ_bands = size(spins,2);
            num_ks = size(spins,3)
            els = zeros(1,num_ks);
            for j = 1:num_ks
                if j ~= num_ks
                    next = j+1;
                else
                    next = 1;
                end
                els(j) = det(spins(:,:,j)' * spins(:,:,next));
            end
            wilson_loop = mod(sum(angle(els))/(2*pi),1);
        end
        
        %******************HAMILTONIAN CONSTRUCTION****************
        
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
        
        function ham = hopping_hamiltonian(hops,potentials,sites,open)
            if numel(hops) == 1
                hops = ones(1,sites).*hops;
            end
            if numel(potentials) == 1
                potentials = ones(1,sites).*potentials;
            end
            hamhop = TopologicalInsulator.off_diagonal_matrix(1,hops(1:sites),1,open);
            hampot = diag(potentials(1:sites));
            ham = hamhop + hamhop' + hampot;
        end
        
        function hk = bloch_hamiltonian(ham,k,cell_size,max_cell_hopping)
            hk = zeros(cell_size);
            for j = 0:max_cell_hopping
                hk = hk + exp(1i*k*j)*ham(j*cell_size + (1:cell_size),1:cell_size);
                if j ~= 0
                    hk = hk + exp(-1i*k*j)*ham(1:cell_size,j*cell_size + (1:cell_size));
                end
            end
        end
        
        %*********************CORRELATION MATRIX********************
        
        function ent = entanglement_entropy_from_correlation_matrix(corrmat, site1, site2)
            redcorrmat = corrmat(site1:site2,site1:site2);
            evals = eig(redcorrmat);
            evals(abs(evals) < 1.e-9 | abs(evals - 1) < 1.e-9) = [];
            H2 = @(x) -x*log(x) - (1-x)*log(1-x);
            ent = sum(arrayfun(H2,evals));
        end
        
        function spec = entanglement_spectrum_from_correlation_matrix(corrmat, site1, site2)
            redcorrmat = corrmat(site1:site2,site1:site2);
            evals = eig(redcorrmat);
            spec = log(1 - evals) - log(evals);
        end
        
        function v = two_particle_expectation(corrmat,creation_vec,annihilation_vec)
            if size(creation_vec,1) ~= 1 || size(annihilation_vec,1) ~= 1
                error('Wavefunctions to be provided as row vectors');
            end
            v = creation_vec * corrmat * annihilation_vec.';
        end
        
        function [mats_out,ks] = k_space_matrix(mat_in,cell_size)
            sites = size(mat_in,1);
            cells = sites/cell_size;
            mats_out = cell(1,cells);
            ks = 2*pi*(0:(cells-1))/cells;
            
            for j = 1:cells
                k = ks(j);
                bloch = exp(1i * k * (0:(cells-1))).';
                bloch_basis_vectors = zeros(sites,cell_size);
                for p = 1:cell_size
                    b = zeros(cell_size,1);
                    b(p) = 1;
                    bloch_basis_vectors(:,p) = kron(bloch,b);
                end
                mats_out{j} = bloch_basis_vectors' * mat_in * bloch_basis_vectors/cells;
            end
        end
        
        function invar = invariant_from_corrmats(mats,ks)
            assert(iscell(mats),'k space matrices should be provided as cell');
            assert(numel(mats) == numel(ks));
            inv_mats = cell(size(mats));
            for j = 1:numel(mats)
                inv_mats{j} = inv(mats{j});
            end
            
            invar = 0;
            for k_ind = 1:numel(ks)
                next = k_ind + 1;
                prev = k_ind - 1;
                if next > numel(ks)
                    next = next - numel(ks);
                end
                if prev < 1
                    prev = prev + numel(ks);
                end
                deriv = (inv_mats{next} - inv_mats{prev});
                invar = invar + trace(mats{k_ind}*deriv);
            end
            
            invar = invar/(2*pi);
        end
        
        
        %********************WANNIER***************************
        
        function xs_sorted = wannier_centres(corrmat,prev_centres,cell_size)
            sites = size(corrmat,1);
            cells = sites/cell_size;
            pos_values = exp(2i*pi*kron(1:(cells),ones(1,cell_size))/cells);
            %pos_values = exp(((0:(sites-1)))*2i*pi/sites);
            proj_pos = corrmat*diag(pos_values)*corrmat;
            centres = eig(proj_pos);
            xs = angle(centres(abs(centres) > 1.e-6))*cells/(2*pi);
            
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
       
        function c_matrix = current_operator_from_hamiltonian(ham,site,max_range)
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
        
        
        %********************PLOTTING*********************
        
        function fig_handle = plot_entanglement_spectrum(specs1,specs2,times)
            assert(size(times,2) == size(specs1,2),...
                'Dimension of times does not agree with spectrum provided');
            pt_to_inch = 72;
            width = 246; height = 170;
            fig_handle = figure('Name', 'Entanglement spectrum',...
                'Units','inches','Position',[4,3,width/pt_to_inch,height/pt_to_inch]);
            
            subplot(2,1,1)
            
            
            
            plot(times,abs(specs1));
            ax1 = gca;
            ylim([0,0.4]);
            set(ax1,'TickLabelInterpreter','latex');
            set(ax1,'FontSize',8);
            set(ax1,'FontName','Times');
            
            set(ax1,'XTickLabel',[]);
            set(ax1,'YTick',[0:0.2:0.4]);
            
            subplot(2,1,2)
            plot(times,abs(specs2));
            ax2 = gca;
            xlab = xlabel('Time $t$','interpreter','latex');
            ylab = ylabel('Entanglement energy $\epsilon_i$','interpreter','latex');
            ylim([0,0.4]);
            set(ax2,'TickLabelInterpreter','latex');
            set(ax2,'FontSize',8);
            set(ax2,'FontName','Times');
            set(ax2,'YTick',[0,0.2]);
            
            leftmargin = 0.105;
            rightmargin = 0.03;
            topmargin = 0.03;
            bottommargin = 0.13;
            
            pos1 = get(ax1,'Position');
            pos2 = get(ax2,'Position');
            pos1(1) = leftmargin;
            pos2(1) = leftmargin;
            pos2(2) = bottommargin;
            figheight = (1 - topmargin - pos2(2))/2;
            pos2(4) = figheight;
            pos1(4) = figheight;
            pos1(2) = pos2(2) + figheight -0.006;
            pos1(3) = 1-pos1(1)-rightmargin;
            pos2(3) = 1-pos2(1)-rightmargin;
            set(ax1,'Position',pos1);
            set(ax2,'Position',pos2);
            set(xlab,'units','normalized');
            set(ylab,'units','normalized');
            posx= get(xlab,'Position');
            posx(2) = -0.18;
            set(xlab,'Position',posx);
            posy= get(ylab,'Position');
            posy(1) = -0.07;
            posy(2) = 1;
            set(ylab,'Position',posy);
            
            
            
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

