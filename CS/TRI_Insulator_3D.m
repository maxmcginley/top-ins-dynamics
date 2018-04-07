classdef TRI_Insulator_3D < Bloch_Hamiltonian
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        so
        mass
    end
    
    methods
        function obj = TRI_Insulator_3D(so,mass)
            obj.so = so;
            obj.mass = mass;
        end
    end
    
    methods (Access = protected)
        function ham = calculate_bloch_hamiltonian(obj,ks)
            pz = [[1,0];[0,-1]];
            px = [[0,1];[1,0]];
            py = [[0,-1i];[1i,0]];
            ham = obj.so*sin(ks(1))*kron(px,px);
            ham = ham + obj.so*sin(ks(2))*kron(px,py);
            ham = ham + obj.so*sin(ks(3))*kron(px,pz);
            
            ham = ham + obj.mass*kron(pz,eye(2));
        end
    end
end

