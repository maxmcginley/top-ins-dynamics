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
            invar = 0;
            nkv1 = numel(kvals{1}); nkv2 = numel(kvals{2}); nkv3 = numel(kvals{3});
            prevw = warning;
            warning('error','MATLAB:logm:nonPosRealEig');
            for ind_1 = 1:nkv1
                for ind_2 = 1:nkv2
                    for ind_3 = 1:nkv3
                        next1 = mod(ind_1,nkv1) + 1; next2 = mod(ind_2,nkv2) + 1; next3 = mod(ind_3,nkv3) + 1;
                        
                        eAx = spinors{ind_1,ind_2,ind_3}' * spinors{next1,ind_2,ind_3};
                        eAy = spinors{ind_1,ind_2,ind_3}' * spinors{ind_1,next2,ind_3};
                        eAz = spinors{ind_1,ind_2,ind_3}' * spinors{ind_1,ind_2,next3};
                        
%                         eAx = eAx / (abs(det(eAx)));
%                         eAy = eAy / (abs(det(eAy)));
%                         eAz = eAz / (abs(det(eAz)));
                        
                        eFxy = (eAx * spinors{next1,ind_2,ind_3}' * spinors{next1,next2,ind_3}) / ...
                             (eAy * spinors{ind_1,next2,ind_3}' * spinors{next1,next2,ind_3});
                        eFyz = (eAy * spinors{ind_1,next2,ind_3}' * spinors{ind_1,next2,next3}) / ...
                             (eAz * spinors{ind_1,next2,ind_3}' * spinors{ind_1,next2,next3});
                        eFzx = (eAz * spinors{ind_1,ind_2,next3}' * spinors{next1,ind_2,next3}) / ...
                             (eAx * spinors{next1,ind_2,ind_3}' * spinors{next1,ind_2,next3});
                        
                        
                        if any(abs([det(eAx),det(eAy),det(eAz)]) < 1.e-5)
                            warning('Singular link variables');
                        end
                        
                        try
                            Ax = logm(eAx); Ay = logm(eAy); Az = logm(eAz);
                        catch ME
                            disp(['Oh no ',ME.identifier]);
                        end
                        
%                         trace(Ax*logm(eFyz)) + trace(Az*logm(eFxy)) + trace(Ay*logm(eFzx))
                        
                        invar = invar + (-1/(8*(pi^2)))*(0.5*(trace(Ax*logm(eFyz)) + ...
                            trace(Az*logm(eFxy)) + trace(Ay*logm(eFzx))) - trace(Ax*(Ay*Az - Az*Ay)));
                        
                    end
                end
            end
            warning(prevw);
        end
        
        %Determinant of U^P := P U P
%         function d = link_determinant(link)
%             es = eigs(
%         end
        
    end
end

