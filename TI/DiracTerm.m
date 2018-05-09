classdef DiracTerm < DiracHamiltonian
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        basis
        name
        k_function
        r_function
    end
    
    methods (Access = protected)
        function obj = DiracTerm(parent,name,k_dependence)
            obj@DiracHamiltonian(parent.dimension_names,parent.sector_names);
            obj.name = name;
            assert(isa(parent,'DiracHamiltonian'));
            obj.parse_k_dependence(k_dependence);
        end
        
        function parse_k_dependence(obj,k_dependence)
            if isempty(k_dependence)
                obj.k_function = @(x,y,z) 1;
                return;
            end
            strs = strsplit(k_dependence);
            assert(mod(numel(strs),2) == 0,'Must have even number of strings');
            obj.k_function = NaN(numel(strs)/2,2);
            for j = 1:(numel(strs)/2)
                assert(strs{2*j}(1) == 'k');
                idx = find(ismember(obj.dimension_names,strs{2*j}(2:end)));
                if isempty(idx)
                    error('Unrecognized dimension name');
                end
                if idx == obj.dimensions
                    if ~isempty(obj.r_function)
                        error('Object`s r_function already set');
                    end
                    switch strs{2*j - 1}
                        case 'cos'
                            obj.r_function = +1;
                        case 'sin'
                            obj.r_function = -1;
                        otherwise
                            error('Unrecognized function of k');
                    end
                else
                    obj.k_function(j,2) = idx;
                    switch strs{2*j - 1}
                        case 'cos'
                            obj.k_function(j,1) = +1;
                        case 'sin'
                            obj.k_function(j,1) = -1;
                        otherwise
                            error('Unrecognized function of k');
                    end
                end
            end
        end
        
        function ham = evaluate(obj,sites,open)
        end
    end
    
end

