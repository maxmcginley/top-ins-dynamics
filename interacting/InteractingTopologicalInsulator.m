classdef InteractingTopologicalInsulator < EDSystem
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimUp
        dimDown
        S
    end
    
    methods
        function f = InteractingTopologicalInsulator(sites,fill,tspin,S)
            if S && (mod(fill+tspin,2) ~= 0 || tspin > fill || fill + tspin > 2*sites) || (~S && tspin ~= fill)
                error('Asking for a spin incompatible with number of fermions');
            end
            
            f.S = S;
            f.sites = sites;
            f.restrs = [fill, tspin];
            f.dimUp = nchoosek(sites,(fill+tspin)/2);
            if S
                f.dimDown = nchoosek(sites,(fill-tspin)/2);
            else
                f.dimDown = 1;
            end
            f.dim = f.dimUp * f.dimDown;
            
            f = f.generate_basis();
        end
        
        
        
        %**************HAMILTONIANS*****************
        
        function [H] = interaction_hamiltonian(obj,J)
            if obj.S 
                int_fun = @(x) J*(x==5);
            
                matels = sum(arrayfun(int_fun,obj.basis),2).';
                rows = 1:obj.dim;

                H = sparse(rows,rows,matels,obj.dim,obj.dim);
            else
                matels = cellfun(@InteractingTopologicalInsulator.spinless_diagonal_term,num2cell(obj.basis,2)).';
                rows = 1:obj.dim;

                H = sparse(rows,rows,matels.*J,obj.dim,obj.dim);
            end
        end
        
        function [H] = kinetic_hamiltonian(obj,Js)
            if numel(Js) == 1
                Js = ones(1,obj,sites).*Js;
            end
            
            matels = zeros(obj.dim*obj.sites*2,1);
            rows = zeros(obj.dim*obj.sites*2,1);
            cols = zeros(obj.dim*obj.sites*2,1);

            pointer = 0;

            for i = 1:obj.dim
                vec_in = obj.basis(i,:);
                for j = 1:obj.sites
                    cr = j;
                    ann = uint32(mod(j,obj.sites) + 1);
                    if (InteractingTopologicalInsulator.cont_spin(vec_in,ann,2)) && ~(InteractingTopologicalInsulator.cont_spin(vec_in,cr,2))
                        [next, m] = obj.gen_vector(vec_in, ann, cr, 2);
                        %disp(['Hopping from ',num2str(ann),' to ', num2str(cr), ' in basis state ', num2str(i), ' to basis state ', num2str(obj.ind(next)), ' matrix element ', num2str(m*t)]);
                        pointer = pointer + 1;

                        matels(pointer,1) = m*Js(j);
                        rows(pointer,1) = i;
                        cols(pointer,1) = obj.ind(next);

                        pointer = pointer + 1;

                        %disp(['Hopping from ',num2str(cr),' to ', num2str(ann), ' in basis state ', num2str(obj.ind(next)), ' to basis state ', num2str(i), ' matrix element ', num2str(m*conj(t))]);
                        
                        matels(pointer,1) = m*conj(Js(j));
                        rows(pointer,1) = obj.ind(next);
                        cols(pointer,1) = i;
                    end

                    if (InteractingTopologicalInsulator.cont_spin(vec_in,ann,3)) && ~(InteractingTopologicalInsulator.cont_spin(vec_in,cr,3))
                        [next, m] = obj.gen_vector(vec_in, ann, cr, 3);
                        pointer = pointer + 1;

                        matels(pointer,1) = m*Js(j);
                        rows(pointer,1) = i;
                        cols(pointer,1) = obj.ind(next);

                        pointer = pointer + 1;

                        matels(pointer,1) = m*conj(Js(j));
                        rows(pointer,1) = obj.ind(next);
                        cols(pointer,1) = i;
                    end
                end
            end

            H = sparse(rows(1:pointer,1),cols(1:pointer,1),matels(1:pointer,1),obj.dim,obj.dim);
        end
        
        %****************STATES*******************
        
        function vec_out = apply_fermion_bilinear_to_state(obj, ann_site, cr_site, spin, vec_in)
            assert(size(vec_in,1) == obj.dim,'Vector in must be of Hilbert space dimension');
            vec_out = zeros(obj.dim,size(vec_in,2));
            for j = 1:obj.dim
                basis_vec_in = obj.basis(j,:);
                if InteractingTopologicalInsulator.bilinear_will_annihilate(basis_vec_in,ann_site,cr_site,spin)
                    continue;
                end
                [next,coeff] = obj.gen_vector(basis_vec_in, ann_site, cr_site, spin);
                vec_out(obj.ind(next),:) = vec_out(obj.ind(next),:) + coeff*vec_in(j,:);
            end
        end
        
        %Returns the matrix of <c^\dagger _j c_k> wrt the state
        function corrmat = correlation_matrix(obj,state)
            if obj.S
                error('Spinful correlation matrix not supported');
            else
                corrmat = zeros(obj.sites);
                for j = 1:obj.sites
                    for k = 1:obj.sites
                        dot_state = obj.apply_fermion_bilinear_to_state(k,j,2,state);
                        corrmat(j,k) = state' * dot_state;
                    end
                end
            end
        end
        
        %***********HASHING METHODS******************
        
        function obj = generate_basis(obj)
            fill = obj.restrs(1);
            spin = obj.restrs(2);
            
            upstart = [ones(1,(fill+spin)/2).*2, zeros(1,obj.sites - (spin+fill)/2)];
            if obj.S
                downstart = [ones(1,(fill-spin)/2).*3, zeros(1,obj.sites - (-spin+fill)/2)];
            else
                downstart = zeros(1,obj.sites);
            end

            obj.basis = zeros(obj.dimUp*obj.dimDown,obj.sites);

            upcurr = upstart;
            i = 1;

            for iUp = 1:obj.dimUp
                downcurr = downstart;
                obj.basis(i,:) = downcurr + upcurr;
                i = i+1;
                for iDown = 1:(obj.dimDown-1)
                    [downcurr, term] = InteractingTopologicalInsulator.advance_sub(downcurr,obj.sites);
                    if term
                        error('Down-spin subbasis reached end before dimDown reached');
                    end
                    obj.basis(i,:) = downcurr + upcurr;
                    i = i+1;
                end
                [upcurr, term] = InteractingTopologicalInsulator.advance_sub(upcurr,obj.sites);
                if term && iUp ~= obj.dimUp
                    error('Up-spin subbasis reached end before dimDown reached');
                end
            end

            [obj.hashes,obj.ind] = sort(cellfun(@EDSystem.basis_hash,num2cell(obj.basis,2)));
        end
        
        
        
        function [ index, coeff ] = gen_vector(obj, vec_in, ann, cr, spin )
            if ~obj.valid_spin(vec_in(1,ann))
                index = -1;
                coeff = nan(1);
                warning('gen_vector has returned nan');
            else
               hash = EDSystem.basis_hash(vec_in);
               hash = hash + spin*sqrt(100*double(cr) +29) - spin*sqrt(100*double(ann) + 29);
               index = EDSystem.bin_search(obj.hashes,hash);

               coeff = obj.fermion_sign(vec_in,ann,cr,spin);
            end
        end
        
        function [ret] = valid_spin(~, x)
            ret = (x == 0 || x == 2 || x == 3 || x == 5);
        end
        
        function [ret] = occupation(~,x)
            if x == 5
                 ret = 2;
             elseif x == 2 || x == 3
                 ret = 1;
             elseif x == 0
                 ret = 0;
             else
                 error('occupation - invalid spin number given');
             end
        end
        
        function [sign] = fermion_sign(obj,vec,ann,cr,spin)
            counter = 0;
            
            

            if cr == ann
                sign = 1;
                return;
            end
            if cr < ann
                counter = counter + InteractingTopologicalInsulator.cont_spin(vec, cr,spin);
                [cr, ann] = deal(ann,cr);
            end

            if spin == 2 && InteractingTopologicalInsulator.cont_spin(vec,ann,3)
                counter = counter + 1;
            end
            if spin == 3 && InteractingTopologicalInsulator.cont_spin(vec,cr,2)
                counter = counter + 1;
            end

            for i = (ann+1):(cr-1)
                counter = counter + obj.occupation(vec(1,i));
            end

            sign = (-1).^counter;
        end
        
        function b = one_site_basis(obj)
            if obj.S
                b = [0 2 3 5];
            else
                b = [0 2];
            end
        end
    end
    
    methods(Static)
        function [out, term] = advance_sub(vec,sites)
           nonz = sites;
           pick = 0;
           while (vec(1,nonz) ~= 0 || vec(1,nonz-1) == 0) 
               if nonz <= 2
                   term = 1;
                   out = [];
                   return;
               end
               if vec(1,nonz) ~= 0
                   pick = pick + 1;
               end
               nonz = nonz - 1;
           end

           out = zeros(1,sites);
           out(1,1:nonz) = vec(1,1:nonz);
           out([nonz-1 nonz]) = vec([nonz nonz-1]);
           out(1,(nonz+1):(nonz+pick)) = ones(1,pick).*vec(1,nonz-1);
           term = 0;

        end
        
        function a = bilinear_will_annihilate(vec_in,ann,cr,spin)
            a = false;
            if ~InteractingTopologicalInsulator.cont_spin(vec_in, ann,spin)
                a = true;
                return;
            end
            if InteractingTopologicalInsulator.cont_spin(vec_in, cr,spin) && cr ~= ann
                a = true;
                return;
            end
        end
        
        function [ret] = cont_spin(vec,ind, spin)
            x = vec(ind);
            if spin == 3
                ret = ((x == 3) || (x ==5));
            elseif spin == 2
                ret = ((x == 2) || (x == 5));
            elseif spin == 5
                ret = (x == 5);
            else
                error('cont_spin - invalid spin number given');
            end
        end
        
        function matel = spinless_diagonal_term(vec)
            matel = 0;
            for i=1:numel(vec)
                plus = i + 2;
                if i == numel(vec)
                    plus = 2;
                elseif i == numel(vec)-1
                    plus = 1;
                end
                matel = matel + (vec(i))*(vec(plus))/4;
            end
        end
        
        function matel = diagonal_term(vec,Js)
            matel = 0;
            inters = numel(vec);
            for i=1:inters
                if vec(i) == 5
                    if numel(Js) == 1
                        matel = matel + Js;
                    else
                        matel = matel + Js(1,i);
                    end
                end
            end
        end
        
    end
    
end

