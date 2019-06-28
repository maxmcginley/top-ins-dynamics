classdef TopologicalInsulator_Ladder < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hoppings
        RF
        raman
        q
        p
        trap
    end
    
    methods
        function obj = TopologicalInsulator_Ladder(sites,hoppings,RF,raman,q,p,trap)
            %q is number of PHYSICAL sites per unit cell
            %Hoppings is a 1 x 2q vector of hopping matrix elements.
            %Neighbouring pairs of elements denote the different hopping
            %strengths on the two legs for the same physical sites
            ham = TopologicalInsulator_Ladder.Ladder_hamiltonian(sites,hoppings,RF,raman,q,p,trap);
            obj@TopologicalInsulator(ham,2*q);
            obj.hoppings = hoppings;
            obj.RF = RF;
            obj.raman = raman;
            obj.q = q;
            obj.p = p;
            obj.trap = trap;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
            intras = zeros(1,2*obj.q - 1);
            intras(1:2:end) = obj.raman*exp(2i*pi*(0:(obj.q-1))*obj.p/obj.q) + obj.RF;
            ham_k = diag(intras,1);
            ham_k = ham_k + diag(obj.hoppings(1:(end-2)),2);
            ham_k(2*obj.q - 1,1) = obj.hoppings(end-1)*exp(1i*k);
            ham_k(2*obj.q,2) = obj.hoppings(end)*exp(1i*k);
            ham_k = ham_k + ham_k';
        end
    end
    
    methods (Access = private)
%         function ham = flux_ham(obj)
%             diags = arrayfun(@(x) TopologicalInsulator_Ladder.BL_diag_entry(x,obj.p/2,gamma,g),0:(N-1));
%             offsup = conj([repmat([obj.hopping1,obj.hopping2],1,(N/2)-1),t1]);
%             offsdown = [repmat([t2,t1],1,(N/2)-1),t2];
%             ham = diag(diags) + diag(offsup,1) + diag(offsdown,-1);
%             ham(1,N) = conj(t1)*exp(-1i*k);
%             ham(N,1) = t2*exp(1i*k);
% 
%             UNCON_CELL = false;
% 
%             if UNCON_CELL
%                 ham(N,1) = t2;
%                 ham(1,1) = ham(1,1)*exp(1i*k);
%                 ham(1,2) = ham(1,2)*exp(-1i*k);
%             end
%         end

        
    end
    
    methods (Static, Access = private) 
        function d = BL_diag_entry(j,pN,gamma,g)
            if mod(j,2) == 0
                d = gamma*exp(-2i*pi*pN*j) + conj(g);
            else
                d = gamma*exp(2i*pi*pN*j) + g;
            end
        end
    end
    
    methods (Static)
        function ham = Ladder_hamiltonian(sites,hoppings,RF,raman,q,p,trap)
            %q is the number of PHYSICAL sites in a unit cell (2q synthetic
            %sites)
            assert(mod(sites,q) == 0, 'Must have commensurate number of sites');
            %Hoppings: a 1 x 2[q] vector of hopping matrix elements between
            %adjascent sites on the same ladder legs
            assert(numel(hoppings) == 2*q,...
                'Number of inter-site hopping matrix elements must match unit cell');
            
            cells = sites/q;
            
            %Hopping within a physical site - raman acquires a phase 2\pi *
            %p / q
            intra_hopping = zeros(1,2*q);
            intra_hopping(1:2:end) = raman*exp(2i*pi*(0:(q-1))*p/q) + RF;
            
            %Harmonic trap
            trap_els = zeros(1,2*sites);
            trap_els(1:2:end) = (((1:sites) - (sites/2)).^2)*trap;
            trap_els(2:2:end) = trap_els(1:2:end);
            
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,intra_hopping,cells,true);
            ham3 = TopologicalInsulator.off_diagonal_matrix(2,hoppings,cells,true);
            ham = ham1 + ham3;
            ham = ham + ham' + diag(trap_els);
        end
        

        
        %*************SYMMETRIES********************
        
        function C = nambu_operator(index,sites)
            if mod(sites,2) ~= 0
                error('Nambu operator not defined on TIs with odd sites');
            end
            switch index
                case 1
                    tau = [[0 1];[1 0]];
                case 2
                    tau = [[0 -1i];[1i 0]];
                case 3
                    tau = [[1 0];[0 -1]];
                otherwise
                    error('Invalid Pauli index');
            end
            C = kron(eye(sites/2),tau);
        end
        
        function chi = test_symmetries(corrmat)
            sites = size(corrmat,1);
            if mod(sites,4) == 0
                chi_op = diag(repmat([1,-1,-1,1],1,sites/4));
            elseif mod(sites,4) == 2
                chi_op = diag([repmat([1,-1,-1,1],1,(sites-2)/4),1,-1]);
            else
                error('Incompatible matrix size');
            end
            
            mat_chi = chi_op * corrmat * chi_op';
            
            chi = sum(sum(abs(corrmat + mat_chi)));
            
        end
    end
end

