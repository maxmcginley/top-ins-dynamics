classdef (Abstract) EDSystem
    %EDSystem Abstract class as a template for boson, fermion, and
    %spin chain systems
    %   Detailed explanation goes here
    
    properties (SetAccess = protected)
        basis
        dim
        sites
        restrs
        hashes
        ind
    end
    
    methods (Abstract)
        generate_basis(obj)
        
        interaction_hamiltonian(obj,J)
        kinetic_hamiltonian(obj,Js)
        
        gen_vector(obj, vec_in, ann, cr, spin )
        valid_spin(obj, x)
        occupation(obj, x)
        
        one_site_basis(obj)
    end
    
    methods (Static)
        function [ hash ] = basis_hash( vec )

            hash = 0;
            for i = 1:numel(vec)
                hash = hash + EDSystem.basis_hash_element(i,vec(1,i));
            end

        end
        
        function h = basis_hash_element(i, val)
            h = sqrt(100*i + 29)*val;
        end
        
        function [v] = gives_non_zero(vec, j, jsp, k, ksp)
            v = (vec(1,j) ~= jsp) && (vec(1,k) ~= ksp);
        end
        
        function [ ind ] = bin_search( hashes, h )
            %UNTITLED6 Summary of this function goes here
            %   Detailed explanation goes here

            tol = 1.e-7;
            offset = 0;

            while 1
                range = size(hashes,1);
                est = idivide(uint32(range)- 1, 2,'floor') + 1;
                if  abs(h - hashes(est)) < tol
                    ind = est + offset;
                    break
                end

                if h > hashes(est)
                    hashes = hashes((est+1):range,1);
                    offset = offset + est;
                else
                    hashes = hashes(1:(est-1),1);
                end

                if h < hashes(1,1) - tol || h > hashes(size(hashes,1),1)+ tol
                    ind = -1;
                    warning('Binary search could not find hash');
                    break;
                end
            end

        end
    end
    
    methods
        
        function val = apply_one_site_operator (obj,op,site,fullvec)
            val = 0;
            osb = obj.one_site_basis;
            for i = 1:obj.dim
                if obj.basis(i,site) ~= osb(1)
                    continue;
                end
                subbasis = zeros(numel(osb),1);
                subbasis(1,1) = i;
                h1 = EDSystem.basis_hash(obj.basis(i,:));
                for j = 2:numel(osb)
                    hash = h1 + EDSystem.basis_hash_element(site,osb(j) - osb(1));
                    subbasis(j,1) = obj.ind(EDSystem.bin_search(obj.hashes,hash));
                end
                subvec = fullvec(subbasis,1);
                val = val + subvec' * op * subvec;
            end
        end
        
        function val = apply_many_site_operator (obj, ops, subsites, fullvec)
            val = 0;
            osb = obj.one_site_basis;
            
            function hash_diff = sub_occ(j)
                occups = zeros(1,numel(subsites));
                hash_diff = 0;
                for k = 1:numel(subsites)
                    occups(1,k) = floor(mod(j-1,numel(osb)^k)/(numel(osb)^(k-1))) + 1;
                    hash_diff = hash_diff + EDSystem.basis_hash_element(subsites(k),osb(occups(1,k)) - osb(1));
                end
                %display(['Index j = ', num2str(j), ' gives occupations ', num2str(occups)]);
            end
            
            function A = full_op(ops)
                if size(ops,3) == 1
                    A = ops(:,:,1);
                else
                    A = kron(full_op(ops(:,:,2:size(ops,3))), ops(:,:,1));
                end
            end
            
            full = full_op(ops);
            
            for i = 1:obj.dim
                if ~isequal(obj.basis(i,subsites), ones(1,numel(subsites))*(osb(1)))
                    continue;
                end
                subbasis = zeros(numel(osb)^numel(subsites),1);
                subbasis(1,1) = i;
                
                h1 = EDSystem.basis_hash(obj.basis(i,:));
                for j = 2:numel(osb)^numel(subsites)
                    hash = h1 + sub_occ(j);
                    subbasis(j,1) = obj.ind(EDSystem.bin_search(obj.hashes,hash));
                end
                subvec = fullvec(subbasis,1);
                val = val + subvec' * full * subvec;
            end
        end
        
    end
    
    
    
end

