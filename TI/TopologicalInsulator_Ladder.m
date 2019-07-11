classdef TopologicalInsulator_Ladder < TopologicalInsulator
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hoppings
        RF
        raman
        q %Periodicity of the flux cell
        p %Flux per unit cell is 2\pi p/q
        N %Number of internal states per site
        trap
    end
    
    methods
        function obj = TopologicalInsulator_Ladder(sites,hoppings,RF,raman,q,p,N,trap)
            %q is number of PHYSICAL sites per unit cell
            %Hoppings is a 1 x (Nq) vector of hopping matrix elements.
            %Neighbouring pairs of elements denote the different hopping
            %strengths on the two legs for the same physical sites
            ham = TopologicalInsulator_Ladder.Ladder_hamiltonian(sites,hoppings,RF,raman,q,p,N,trap);
            obj@TopologicalInsulator(ham,N*q);
            obj.hoppings = hoppings;
            obj.RF = RF;
            obj.raman = raman;
            obj.q = q;
            obj.p = p;
            obj.N = N;
            obj.trap = trap;
        end
        
        
        function ham_k = BL_k_hamiltonian(obj,k)
            [intras,inters] = TopologicalInsulator_Ladder.get_hopping_elements(...
                obj.hoppings,obj.RF,obj.raman,obj.q,obj.p,obj.N);
            ham_k = diag(intras(1:end-1),1);
            ham_k = ham_k + diag(inters(1:(end-obj.N)),obj.N);
            ham_k((obj.N*(obj.q - 1) + 1):end,1:obj.N) = ...
               ham_k((obj.N*(obj.q - 1) + 1):end,1:obj.N) + ...
               diag(inters((end-obj.N+1):end))*exp(1i*k);
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
        
        function [intras,inters] = get_hopping_elements(hoppings,RF,raman,q,p,N)
            line = [reshape(raman,N-1,1);0];
            phases = exp(2i*pi*(0:(q-1))*p/q);
            intras = line * phases;
            intras = intras + ([ones(N - 1,1)*RF;0] * ones(1,q));
            intras = reshape(intras,1,q*N);
            inters = reshape(ones(N,1)*reshape(hoppings,1,q),1,q*N);
        end
    end
    
    methods (Static)
        function ham = Ladder_hamiltonian(sites,hoppings,RF,raman,q,p,N,trap)
            %q is the number of PHYSICAL sites in a unit cell (2q synthetic
            %sites)
            assert(mod(sites,q) == 0, 'Must have commensurate number of sites');
            %Hoppings: a 1 x 2[q] vector of hopping matrix elements between
            %adjascent sites on the same ladder legs
            assert(numel(hoppings) == q,...
                'Number of inter-site hopping matrix elements must match unit cell');
            assert(numel(raman) == N-1,'Must specify all Raman couplings');
            cells = sites/q;
            
            [intras,inters] = TopologicalInsulator_Ladder.get_hopping_elements(...
                hoppings,RF,raman,q,p,N);
            
            %Hopping within a physical site - raman acquires a phase 2\pi *
            %p / q
            
            %Harmonic trap
            trap_els = zeros(1,N*sites);
            tr = (((1:sites) - ((sites+1)/2)).^2)*trap;
            trap_els = kron(tr,ones(1,N));
            
            ham1 = TopologicalInsulator.off_diagonal_matrix(1,intras,cells,true);
            ham3 = TopologicalInsulator.off_diagonal_matrix(N,inters,cells,true);
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
        
        function chi = test_symmetries(corrmat,N)
            sites = size(corrmat,1)/N;
            line = ones(1,N);
            line(2:2:end) = -1;
            if mod(sites,2) == 0
                chi_op = diag(repmat([line,-line],1,sites/2));
            else
                chi_op = diag([repmat([line,-line],1,(sites-1)/2),line]);
            end
            
            mat_chi = chi_op * corrmat * chi_op';
            
            chi = sum(sum(abs(corrmat + mat_chi)));
            
        end
    end
end

