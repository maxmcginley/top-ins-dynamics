classdef DiracHamiltonian
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bands
        dimension
        dimension_names
        sector_names
    end
    
    properties (GetAccess = protected, SetAccess = protected)
        terms
    end
    
    methods
        function obj = DiracHamiltonian(dimension_names,sector_names)
            [obj.bands,obj.dimension] = DiracHamiltonian.parse_names(dimension_names,sector_names);
            obj.terms = cell();
            obj.dimension_names = dimension_names;
            obj.sector_names = sector_names;
        end
        
        function ham = hamiltonian(obj,sites,open)
            ham = zeros(sites);
            for j = 1:numel(obj.terms)
                ham = ham + obj.terms{j}.evaluate(sites,open);
            end
        end
    end
    
    methods (Static,Access = protected)
        function [bands,dimensions] = parse_names(dimension_names,sector_names)
            assert(iscell(sector_names),'Sector names to be provided as cell');
            assert(iscell(dimension_names),'Dimension names to be provided as cell');
            bands = 2.^(numel(sector_names));
            dimension = numel(dimension_names);
        end
    end
end

