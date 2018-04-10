classdef (Abstract) Bloch_Hamiltonian
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimension
        bands
    end
    
    methods (Abstract,Access = protected)
        
        ham = calculate_bloch_hamiltonian(obj,ks)
        
    end
    
    methods
        function obj = Bloch_Hamiltonian(dimension,bands)
            obj.dimension = dimension;
            obj.bands = bands;
        end
        
        function ham = hamiltonian(obj,ks)
            assert(numel(ks) == obj.dimension);
            ham = obj.calculate_bloch_hamiltonian(ks);
            assert(size(ham,1) == size(ham,2));
            assert(size(ham,1) == obj.bands);
        end
        
        function [spinors,ens] = ground_state_spinors(obj,ks)
            assert(numel(ks) == obj.dimension);
            k_dims = cellfun(@numel,ks);
            spinors = cell(k_dims);
            num_kpoints = prod(k_dims);
            
            num_bands = NaN;
            
            ens = cell(k_dims);
            
            for k_index = 1:num_kpoints
                [m1,m2,m3] = ind2sub(k_dims,k_index);
                k = [ks{1}(m1),ks{2}(m2),ks{3}(m3)];
                ham_k = obj.hamiltonian(k);
                [evec,d] = eig(ham_k);
                eval = real(diag(d)); occs = eval < 0.0;
                if any(occs ~= [1;1;0;0])
                    error('OCCS');
                end
                if isnan(num_bands)
                    num_bands = sum(occs);
                else
                    if num_bands ~= sum(occs)
                        error('Inconsistent number of bands');
                    end
                end
                if any(abs(eval) < 1.e-4)
                    error('Found singular Hamiltonian');
                end
                spinors{m1,m2,m3} = evec(:,occs);
                ens{m1,m2,m3} = eval(occs);
            end
        end
    end
    
    methods (Static)
        
        function ks = generate_kvals(k_step,dim)
            kvals = [(-0.5+k_step):k_step:0.5]*(2*pi);
            ks = cell(1,dim);
            for j = 1:dim
                ks{j} = kvals;
            end
        end
        
    end
end

