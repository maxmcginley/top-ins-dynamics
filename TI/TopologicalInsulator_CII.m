classdef TopologicalInsulator_CII < TopologicalInsulator
    %Class CII insulator in 3 dimensions
    %   Detailed explanation goes here
    
    properties
        velocity
        mass
        inversion
        chi_breaking
        z_cells
    end
    
    methods
        function obj = TopologicalInsulator_CII(velocity,mass,inversion,chi_breaking,kx,ky,sites,open)
            ham = TopologicalInsulator_CII.CII_hamiltonian(velocity,mass,inversion,chi_breaking,kx,ky,sites,open);
            obj@TopologicalInsulator(ham,4);
            obj.velocity = velocity; obj.mass = mass; obj.inversion = inversion;
            obj.z_cells = sites; obj.chi_breaking = chi_breaking;
        end
        
        function ham_k = BL_k_hamiltonian(~,~)
            error('Not supported');
        end
        
        
    end
    
    methods(Static)
        %function CII_Hamiltonian: constructs a CII topological insulator
        %Hamiltonian which is fourier transformed in x and y directions
        %only. The ordering of indices is first cell number (z), then tau
        %index (1/2), then spin index (u/d). i.e. [1u,1d,2u,2d, 1u(z+1),
        %...] Note: sites here refers to the number of cells in the z
        %direction
        function ham = CII_hamiltonian(velocity,mass,inversion,chi_breaking,kx,ky,sites,open)
            %mass_new = inversion*(1- 0.25*cos(kx) - 0.25*cos(ky));
            
            mass_new = mass - inversion*(cos(kx) + cos(ky));
            
            %hoppings1 = [1u-1d, 1d-2u, 2u-2d, 2d-1u(+1)];
            hoppings_1 = [velocity*(sin(kx)-1i*sin(ky)),0,-velocity*(sin(kx)-1i*sin(ky)),0];
            
            %hoppings2 = [1u-2u, 1d-2d, 2u-1u(+1), 2d-1d(+1)];
            hoppings_2 = [-1i*(mass_new) + chi_breaking*sin(kx),-1i*(mass_new) + chi_breaking*sin(kx),1i*inversion,1i*inversion];
            
            %hoppings3 = [1u-2d, 1d-1u(+1), 2u-1d(+1), 2d-2u(+1)];
            %hoppings_3 = inversion*[0,1i,0,-1i];
            
            %hoppings4 = [1u-1u,1d-1d,2u-2u,2d-2d] all (+1);
            hoppings_4 = velocity*[1i,-1i,-1i,1i];
            
            %hoppings5 = [1u-1d(+1), 1d-2u(+1), 2u-2d(+1), 2d-1u(+2)];
            %hoppings_5 = inversion*[-1i,0,1i,0];
            
            %hoppings6 = [1u-2u(+1), 1d-2d(+1), 2u-1u(+2), 2d-1d(+2)];
            hoppings_6 = [-1i*inversion, -1i*inversion,0,0];
            
            
            
            hop1 = TopologicalInsulator.off_diagonal_matrix(1,hoppings_1,sites,open);
            hop2 = TopologicalInsulator.off_diagonal_matrix(2,hoppings_2,sites,open);
            %hop3 = TopologicalInsulator.off_diagonal_matrix(3,hoppings_3,sites,open);
            hop4 = TopologicalInsulator.off_diagonal_matrix(4,hoppings_4,sites,open);
            %hop5 = TopologicalInsulator.off_diagonal_matrix(5,hoppings_5,sites,open);
            hop6 = TopologicalInsulator.off_diagonal_matrix(6,hoppings_6,sites,open);
            
            ham = hop1 + hop2 + hop4 + hop6;
            ham = ham + ham';
        end
        
        function chi = test_chiral_symmetry(matrix)
            sites = size(matrix,1); cells = sites/4;
            chi_op = kron(eye(cells),kron([[0,1];[1,0]],eye(2)));
            mat_chi = chi_op * matrix * chi_op';
            chi = sum(sum(abs(matrix + mat_chi)));
            
        end
        
        
    end
end

