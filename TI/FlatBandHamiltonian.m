classdef FlatBandHamiltonian
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Access = private)
        function obj = FlatBandHamiltonian()
        end
    end
    
    methods(Static)
        %Artificially creates a flat-band Hamiltonian
        function ham = four_cell_hamiltonian(params_in,breaking,chiral_breaking,cells,open)
            if numel(params_in) == 4
                [a,b,c,s] = params_in{:};
                [t,r,u] = FlatBandHamiltonian.dependent_parameters(a,b,c,s,breaking);
            elseif numel(params_in) == 7
                [a,b,c,t,r,s,u] = params_in{:};
                r = r - breaking;
            else
                error('Invalid number of parameters');
            end
            
            hop1 = TopologicalInsulator.off_diagonal_matrix(1,[t,s,u,b'],cells,open);
            hop3 = TopologicalInsulator.off_diagonal_matrix(3,[r,a',c',0],cells,open);
            hop2 = TopologicalInsulator.off_diagonal_matrix(2,ones(1,4)*chiral_breaking,cells,open);
            
            ham = hop1 + hop3 + hop1' + hop3' + hop2 + hop2';
            
        end
        
        function [hc,hs] = cos_sin_vectors(h0,hc_init)
            hc = hc_init - (hc_init.' * h0) * h0/norm(h0)^2;
            hs = cross(h0,hc);
            hs = hs * norm(hc)/norm(hs);
        end
        
        function ham = two_cell_hamiltonian(h0,hc,hs,cells,open)
            c = 0.5*(hc(3) - 1i*hs(3));
            a = 0.5*(hc(1) + hs(2)) + 0.5i*(hc(2) - hs(1));
            b = 0.5*(hc(1) - hs(2)) - 0.5i*(hc(2) + hs(1));
            
%             k_vals = (0:(cells-1))*(2*pi)/cells;
%             es = zeros(2,cells);
%             for j = 1:numel(k_vals)
%                 k = k_vals(j);
%                 es(:,j) = eig([[h0(3)+ c*exp(1i*k) + c'*exp(-1i*k),...
%                     h0(1)-1i*h0(2)+b*exp(1i*k)+a'*exp(-1i*k)];...
%                     [h0(1)+1i*h0(2)+b'*exp(-1i*k)+a*exp(1i*k),...
%                     -h0(3) - c*exp(1i*k) - c'*exp(-1i*k)]]);
%             end
%             es
            
            hop0 = TopologicalInsulator.off_diagonal_matrix(0,[h0(3),-h0(3)]/2,cells,open);
            hop1 = TopologicalInsulator.off_diagonal_matrix(1,[h0(1)-1i*h0(2),a],cells,open);
            hop2 = TopologicalInsulator.off_diagonal_matrix(2,[c,-c],cells,open);
            hop3 = TopologicalInsulator.off_diagonal_matrix(3,[b,0],cells,open);
            
            ham = hop0 + hop1 + hop2 + hop3;
            ham = ham + ham';
            eig(ham(1:2,1:2) + ham(3:4,1:2) + ham(1:2,3:4))
        end
        
        function [t,r,u] = dependent_parameters(a,b,c,s,breaking)
            t = (c')*a*(b*b' - s*s')/(s'*(a*a' + b*b'));
            r = -(c')*b*(a*a' + s*s')/(s'*(a*a' + b*b')) -breaking;
            u = b*s'/a ;
        end
        
        function [en1, en2] = energies_from_parameters(params_in)
            if numel(params_in) == 4
                [a,b,c,s] = params_in{:};
                [t,r,u] = FlatBandHamiltonian.dependent_parameters(a,b,c,s,0);
            elseif numel(params_in) == 7
                [a,b,c,t,r,s,u] = params_in{:};
            else
                error('Invalid number of parameters');
            end
            sumsq = (a+t)*(a+t)' + (b+r)*(b+r)' + (c+s)*(c+s)' + u*u';
            miscterm = -4*abs((a+t)*u - (b+r)*(c+s)').^2;
            en1 = sqrt(sumsq - sqrt(sumsq^2 + miscterm));
            en2 = sqrt(sumsq + sqrt(sumsq^2 + miscterm));
        end
        
        function params_out = commensurate_parameters(en1,en2,a,b,c,s)
            en1 = abs(en1); en2 = abs(en2);
            if en1 > en2
                [en1,en2] = deal(en2,en1);
            end
            ratio = en2/en1;
            mu = a*a' + b*b' + s*s' + b*b'*s*s'/(a*a');
            beta = abs(b)*(a*a' + s*s')/(abs(s)*(a*a' + b*b'));
            alpha = 1 + (beta)^2 + (a*a'/(s*s'))'*(((b*b' - s*s')/(b*b' + a*a'))^2);
            max_rat = (alpha + sqrt(alpha^2-beta^2))/(beta);
            min_rat = (alpha - sqrt(alpha^2-beta^2))/(beta);
%             if ratio > min_rat && ratio < max_rat
%                 fprintf('Ratio must not be in [%d,%d]',min_rat,max_rat);
%                 error('Cannot achieve desired energy ratio with given parameters');
%             end
%             csq = 2*mu*ratio*(beta*(1 + ratio^2) - 2*alpha*ratio)/...
%                 (beta*(1+ratio^2)^2 + 4*alpha^2*ratio^2);
            csq = ratio*mu/(-alpha*ratio + beta*(1 + ratio^2));
            c = exp(1i*angle(c))*sqrt(csq);
            [t,r,u] = FlatBandHamiltonian.dependent_parameters(a,b,c,s,0);
            [en1_t, en2_t] = FlatBandHamiltonian.energies_from_parameters({a,b,c,t,r,s,u});
            
            rescale_1 = en1/en1_t;
            rescale_2 = en2/en2_t;
            if abs(rescale_1 - rescale_2) > 1.e-9
                error('Did not find energies in desired ratio');
            end
            params_out = num2cell(sqrt(2)*rescale_1*[a,b,c,t,r,s,u]);
        end
    end
end

