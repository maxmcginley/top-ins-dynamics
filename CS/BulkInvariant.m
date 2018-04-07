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
            assert(all(size(spinors) == size(kvals)));
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
                        link1 = det(spinors{ind_1,ind_2,ind_3}' * spinors{next1,ind_2,ind_3}); %U^P_1
                        link2 = det(spinors{ind_1,ind_2,ind_3}' * spinors{ind_1,next2,ind_3}); %U^P_2
                        link3 = det(spinors{ind_1,ind_2,ind_3}' * spinors{ind_1,ind_2,next3}); %U^P_3
                        
                        link4 = det(spinors{next1,ind_2,ind_3}' * spinors{next1,next2,ind_3});
                        link5 = det(spinors{ind_1,next2,ind_3}' * spinors{next1,next2,ind_3});
                        
                        exp_curv = (link1*link4)/(link2*link5);
                        
                        invar = invar + angle(exp_curv)*angle(link3) + (2/3)*angle(link1)*angle(link2)*angle(link3);
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

