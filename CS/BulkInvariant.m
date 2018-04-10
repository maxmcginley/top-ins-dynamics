classdef BulkInvariant
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dimension
        kvals
        spinors
        bands
    end
    
    methods
        function obj = BulkInvariant(dimension,kvals,spinors)
            obj.dimension = dimension;
            assert(iscell(kvals)); assert(numel(kvals) == dimension);
            obj.kvals = kvals;
            obj.dimension = dimension;
            assert(iscell(spinors));
            assert(all(size(spinors) == cellfun(@numel,kvals)));
            obj.spinors = spinors;
            obj.bands = size(spinors{1},1);
        end
        
        function invar = invariant(obj)
            switch obj.dimension
                case 3
                    invar = BulkInvariant.ThreeD_CS(obj.kvals,obj.spinors);
                otherwise
                    error('Dimension not supported');
            end
        end
    end
    
    methods (Static)
        
        function invar = ThreeD_CS(kvals,spinors)
            perm1 = [2 3 1];
            perm2 = [3 1 2];
            
            invar = 2*BulkInvariant.ThreeD_CS_component(kvals,spinors);
            invar = invar + 2*BulkInvariant.ThreeD_CS_component(...
                permute(kvals,perm1),permute(spinors,perm1));
            invar = invar + 2*BulkInvariant.ThreeD_CS_component(...
                permute(kvals,perm2),permute(spinors,perm2));
        end
        
        function invar = ThreeD_CS_component(kvals,spinors)
            invar = 0;
            nkv1 = numel(kvals{1}); nkv2 = numel(kvals{2}); nkv3 = numel(kvals{3});
            for ind_1 = 1:nkv1
                for ind_2 = 1:nkv2
                    for ind_3 = 1:nkv3
                        next1 = mod(ind_1,nkv1) + 1; next2 = mod(ind_2,nkv2) + 1; next3 = mod(ind_3,nkv3) + 1; 
                        link(1) = det(spinors{ind_1,ind_2,ind_3}' * spinors{next1,ind_2,ind_3}); %U^P_1
                        link(2) = det(spinors{ind_1,ind_2,ind_3}' * spinors{ind_1,next2,ind_3}); %U^P_2
                        link(3) = det(spinors{ind_1,ind_2,ind_3}' * spinors{ind_1,ind_2,next3}); %U^P_3
                        
                        link(4) = det(spinors{next1,ind_2,ind_3}' * spinors{next1,next2,ind_3});
                        link(5) = det(spinors{ind_1,next2,ind_3}' * spinors{next1,next2,ind_3});
                        
                        if any(abs(link) < 1.e-5)
                            warning('Singular link variables');
                        end
                        
                        exp_curv = (link(1)*link(4))/(link(2)*link(5));
                        
                        invar = invar + (-1/(8*(pi^2)))*(angle(exp_curv)*angle(link(3)) + (2/3)*angle(link(1))*angle(link(2))*angle(link(3)));
                    end
                end
            end
        end
        
        %Determinant of U^P := P U P
%         function d = link_determinant(link)
%             es = eigs(
%         end
        
    end
end

