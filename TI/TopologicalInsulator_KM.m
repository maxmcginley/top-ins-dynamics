classdef TopologicalInsulator_KM < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hopping
        so
        y_cells
    end
    
    methods
        function obj = TopologicalInsulator_KM(hopping,so,rashba,sublattice,kx,sites,open)
            ham = TopologicalInsulator_KM.KM_hamiltonian(hopping,so,rashba,sublattice,kx,sites,open);
            obj@TopologicalInsulator(ham,8);
            obj.hopping = hopping; obj.so = so;
            obj.y_cells = sites;
        end
        
        function ham_k = BL_k_hamiltonian(~,~)
            error('Not supported');
        end
        
        
    end
    
    methods(Static)
        %function KM_Hamiltonian: constructs a Kane-Mele Hamiltonian which
        %is fourier transformed in one direction only (along the point of
        %the hexagons). Spin labels are carried in the indices. The
        %ordering of indices is first cell number (y), then site index
        %(a-d), then spin (^ or v). i.e.
        %[a^, av, b^, bv, c^, cv, d^, dv, a(y+1)^, ...]
        %Note: sites here refers to the number of cells in the y direction
        function ham = KM_hamiltonian(hopping,so,r,sublattice,k,sites,open)
            %Notation: ^ = up spin, v = down spin. Letters a-d label the
            %4 sites within a cell.
            ek = exp(1i*k);
            %hoppings2 = [a^-b^, av-bv, b^-c^, bv-cv, c^-d^, cv-dv, ...
            %..., d^-a^(+1), dv-av (+1)];
            hoppings_2 = hopping*[1,1,1,1,1,1,exp(1i*k),exp(1i*k)];
            
            %hoppings4 = [a^-c^, av-cv, b^-d^, bv-dv, c^-a^(+1), cv-av(+1), ...
            %..., d^-b^(+1), dv-bv (+1)];
            hoppings_4 = 1i*so*[+1+ek',-1-ek',-1-ek',+1+ek',+1+ek,-1-ek,-1-ek,+1+ek];
            
            %hoppings6 = [a^-d^, av-dv, b^-a^(+1), bv-av(+1), c^-b^(+1), cv-bv(+1), ...
            %..., d^-c^(+1), dv-cv (+1)];
            hoppings_6 = hopping*[ek',ek',0,0,1,1,0,0];
            
            %hoppings8 = [a^-a^, av-av, b^-b^, bv-bv, c^-c^, cv-cv, ...
            %..., d^-d^, dv-dv] all (+1);
            hoppings_8 = 1i*so*[-1,+1,+1,-1,-1,+1,+1,-1];
            
            r3 = sqrt(3);
            
            %rahsba1 = [a^-av, av-b^, b^-bv, bv-c^, c^-cv, cv-d^, ...
            %..., d^-dv, dv-a^ (+1)];
            rashba_1 = 1i*r*[0,1i,0,r3-1i,0,1i,0,(r3-1i)*ek];
            
            %rahsba3 = [a^-bv, av-c^, b^-cv, bv-d^, c^-dv, cv-a^(+1), ...
            %..., d^-av(+1), dv-b^ (+1)];
            rashba_3 = 1i*r*[-1i,0,r3+1i,0,-1i,0,(r3+1i)*ek,0];
            
            %rahsba5 = [a^-cv, av-d^, b^-dv, bv-a^, c^-av, cv-b^(+1), ...
            %..., d^-bv(+1), dv-c^ (+1)];
            rashba_5 = 1i*r*[0,(r3+1i)*ek',0,0,0,r3+1i,0,0];
            
            %rahsba7 = [a^-dv, av-a^(+1), b^-av(+1), bv-b^(+1), c^-bv(+1), cv-c^(+1), ...
            %..., d^-cv(+1), dv-d^ (+1)];
            rashba_7 = 1i*r*[(r3-1i)*ek',0,0,0,r3-1i,0,0,0];
            
            on_site_pots = sublattice*[1,1,-1,-1,1,1,-1,-1];
            
            onsite = TopologicalInsulator.off_diagonal_matrix(0,on_site_pots,sites,open);
            hop2 = TopologicalInsulator.off_diagonal_matrix(2,hoppings_2,sites,open);
            hop4 = TopologicalInsulator.off_diagonal_matrix(4,hoppings_4,sites,open);
            hop6 = TopologicalInsulator.off_diagonal_matrix(6,hoppings_6,sites,open);
            hop8 = TopologicalInsulator.off_diagonal_matrix(8,hoppings_8,sites,open);
            
            hop1 = TopologicalInsulator.off_diagonal_matrix(1,rashba_1,sites,open);
            hop3 = TopologicalInsulator.off_diagonal_matrix(3,rashba_3,sites,open);
            hop5 = TopologicalInsulator.off_diagonal_matrix(5,rashba_5,sites,open);
            hop7 = TopologicalInsulator.off_diagonal_matrix(7,rashba_7,sites,open);
            
            
            ham = hop1 + hop2 + hop3 + hop4 + hop5 + hop6 + hop7 + hop8 + onsite/2;
            ham = ham + ham';
        end
        
        function trs = time_reversal_operator(y_cells)
            spin_exch = [[0, 1];[-1,0]];
            trs = kron(eye(y_cells*4),spin_exch);
        end
        
        function t = matrix_has_trs(mat_k,mat_km)
            assert(all(size(mat_k) == size(mat_km)),'Matrices must have same size');
            trs = TopologicalInsulator_KM.time_reversal_operator(size(mat_k,1)/8);
            t_conj = -trs * conj(mat_k) * trs;
            if sum(sum(abs(t_conj - mat_km))) < 1.e-6
                t = true;
            else
                disp(t_conj - mat_km);
                t = false;
            end
        end
        
    end
end

