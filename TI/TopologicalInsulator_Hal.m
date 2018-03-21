classdef TopologicalInsulator_Hal < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hopping
        nn
        y_cells
    end
    
    methods
        function obj = TopologicalInsulator_Hal(hopping,nn,phase,sublattice,kx,sites,open)
            ham = TopologicalInsulator_Hal.Haldane_hamiltonian(hopping,nn,phase,sublattice,kx,sites,open);
            obj@TopologicalInsulator(ham,4);
            obj.hopping = hopping; obj.nn = nn;
            obj.y_cells = sites;
        end
        
        function ham_k = BL_k_hamiltonian(~,~)
            error('Not supported');
        end
        
        
    end
    
    methods(Static)
        %function Haldane_Hamiltonian: constructs a Haldane Hamiltonian
        %which is fourier transformed in one direction only (along the
        %point of the hexagons). The ordering of indices is first cell
        %number (y), then site index (a-d). i.e. [a, b c, d, a(y+1), ...]
        %Note: sites here refers to the number of cells in the y direction
        function ham = Haldane_hamiltonian(hopping,nn,phase,sublattice,k,sites,open)
            %Notation: Letters a-d label the 4 sites within a cell.
            ek = exp(1i*k);
            ph = exp(1i*phase);
            %hoppings1 = [a-b, b-c, c-d, d-a(+1)];
            hoppings_1 = hopping*[1,1,1,exp(1i*k)];
            
            %hoppings2 = [a-c, b-d, c-a(+1), d-b(+1)];
            hoppings_2 = nn*[ph*(1+ek'),ph'*(1+ek'),ph*(1+ek),ph'*(1+ek)];
            
            %hoppings3 = [a-d, b-a(+1), c-b(+1), d-c(+1)];
            hoppings_3 = hopping*[ek',0,1,0];
            
            %hoppings4 = [a-a, b-b, c-c, d-d] all (+1);
            hoppings_4 = nn*[ph',ph,ph',ph];
            
            
            on_site_pots = sublattice*[1,-1,1,-1];
            
            onsite = TopologicalInsulator.off_diagonal_matrix(0,on_site_pots,sites,open);
            hop1 = TopologicalInsulator.off_diagonal_matrix(1,hoppings_1,sites,open);
            hop2 = TopologicalInsulator.off_diagonal_matrix(2,hoppings_2,sites,open);
            hop3 = TopologicalInsulator.off_diagonal_matrix(3,hoppings_3,sites,open);
            hop4 = TopologicalInsulator.off_diagonal_matrix(4,hoppings_4,sites,open);
            
            
            ham = hop1 + hop2 + hop3 + hop4 + onsite/2;
            ham = ham + ham';
        end
        
        
    end
end

