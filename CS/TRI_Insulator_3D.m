classdef TRI_Insulator_3D < Bloch_Hamiltonian
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        so
        mass
        chi
        hop
    end
    
    methods
        function obj = TRI_Insulator_3D(so,mass,hop,chi)
            obj@Bloch_Hamiltonian(3,4);
            obj.so = so;
            obj.mass = mass;
            obj.hop = hop;
            obj.chi = chi;
        end
    end
    
    methods (Access = protected)
        function ham = calculate_bloch_hamiltonian(obj,ks)
            pz = [[1,0];[0,-1]];
            px = [[0,1];[1,0]];
            py = [[0,-1i];[1i,0]];
            
            ham = obj.so*sin(ks(1))*kron(pz,px);
            ham = ham + obj.so*sin(ks(2))*kron(eye(2),py);
            ham = ham + obj.chi*sin(ks(3))*kron(px,px);
            
            ham = ham + (obj.mass - obj.hop*(3 - sum(cos(ks))))*kron(eye(2),pz);
            
        end
    end
end

