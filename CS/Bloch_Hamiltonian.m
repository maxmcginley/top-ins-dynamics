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
        
        function [spinors,bands] = ground_state_spinors(obj,ks)
            assert(numel(ks) == obj.dimension);
            spinors = cell(size(ks));
            num_kpoints = prod(cell2mat(cellfun(numel,ks)));
            
            bands = NaN;
            
            for k_index = 1:num_kpoints
                multi_index = ind2sub(size(ks),k_index);
                k = [ks{1}(multi_index(1)),ks{2}(multi_index(2)),ks{3}(multi_index(3))];
                ham_k = obj.hamiltonian(k);
                [evec,d] = eig(ham_k);
                eval = real(diag(d)); occs = eval < 0.0;
                if isnan(bands)
                    bands = sum(occs);
                else
                    if bands ~= sum(occs)
                        error('Inconsistent number of bands');
                    end
                end
                spinors{multi_index} = evec(:,occs);
            end
        end
    end
end

