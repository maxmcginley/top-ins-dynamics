classdef TimeEvolution_Noise < TimeEvolution
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hams
        fourier_amplitudes
        frequencies
        num_channels
        static_ham
        static_eval
        static_evec
        static_basis
        signal
        times
    end
    
    methods
        function obj = TimeEvolution_Noise(timestep,num_steps,max_exp,hams,spectra,rand_phase,static_ham,static_basis)
            obj@TimeEvolution(timestep,num_steps,max_exp);
            assert(iscell(hams)); assert(iscell(spectra));
            assert(numel(hams) == numel(spectra),...
                'Must provide equal numbers of Hamiltonians and spectra');
            obj.num_channels = numel(hams);
            obj.static_ham = static_ham;
            [obj.static_evec,eval] = eig(static_ham);
            obj.static_eval = diag(eval);
            obj.static_basis = static_basis;
            if static_basis
                obj.hams = cell(size(hams));
                for j = 1:numel(hams)
                    obj.hams{j} = obj.static_evec' * hams{j} * obj.static_evec;
                end
            else
                obj.hams = hams;
            end
            fqs = TimeEvolution_Noise.generate_frequencies(obj.timestep,obj.num_steps);
            [obj.frequencies, obj.fourier_amplitudes] = TimeEvolution_Noise.randomize_amplitudes(...
                    fqs,spectra,rand_phase);
        end
        
        function unit = calculate_static_unitary(obj,time)
            ph_diag = obj.calculate_static_phases(time);
            if obj.static_basis
                unit = diag(ph_diag);
            else
                unit = obj.static_evec * bsxfun(@times,obj.static_evec',ph_diag);
            end
        end
        
        function ph = calculate_static_phases(obj,time)
            ph = exp(-1i * time * obj.static_eval);
        end
        
        function ham = calculate_hamiltonian(obj,t_step,interaction_pic)
%             if isempty(real_number)
%                 error('Realization number not given');
%             end
            ham = zeros(size(obj.static_ham));
            
            for i = 1:obj.num_channels
                %signal_i = imag(obj.signal(i,mod(t_step,size(obj.signal,2))+1));
                signal_i = sum(imag(exp(1i * double(t_step) * obj.timestep * obj.frequencies) .* obj.fourier_amplitudes(i,:)));
                ham = ham + signal_i*obj.hams{i};
            end
            if obj.static_basis
                if numel(interaction_pic) == 1
                    ph = obj.calculate_static_phases(double(interaction_pic)*obj.timestep);
                    ham = (ph * ph') .* ham;
                else
                    ham = ham + diag(obj.static_eval);
                end
            else
                if numel(interaction_pic) == 1
                    unit = obj.calculate_static_unitary(double(interaction_pic)*obj.timestep);
                    ham = unit * ham * unit';
                else
                    ham = ham + obj.static_ham;
                end
            end
        end
        
        function sigs = calculate_signals(obj,t_step,channel_number)
            sigs = imag(obj.signal(channel_number,t_step));
        end
        
        function obj = allocate_phases(obj,times)
            obj.signal = ifft(obj.fourier_amplitudes,[],2)*size(obj.fourier_amplitudes,2);
            obj.times = times;
        end
        
    end
    
    methods (Static)
        
        function amp = cutoff_function(x,cut)
            amp = exp(-(1 + (x.^2 / (cut^2)).^(-4)));
        end
        
        %NORMALISE TO f(0) = 1
        function amp = spec_function(x,wid)
            %amp = (2*wid/pi)*((wid.^2 + x.^2).^(-1));
            amp =  exp(-2*(wid^(-1))*abs(x));
        end
        
        function spectra = generate_poissonians(widths,amplitudes,varargin)
            assert(numel(widths) == numel(amplitudes),...
                'Must provide equal numbers of widths and amplitudes');
            spectra = cell(size(widths));
            if numel(varargin) >= 1
                cutoffs = varargin{1};
                assert(numel(widths) == numel(cutoffs) || isempty(cutoffs));
            else
                cutoffs = [];
            end
            for i = 1:numel(widths)
                if isempty(cutoffs)
                    spectra{i} = @(x) (amplitudes(i).^2).*TimeEvolution_Noise.spec_function(x,widths(i));
                else
                    spectra{i} = @(x) (amplitudes(i).^2).*TimeEvolution_Noise.cutoff_function(x,cutoffs(i)) ... 
                        .* TimeEvolution_Noise.spec_function(x,widths(i));
                end
            end
        end
        
        %Angular frequencies (times 2 pi)
        function fs = generate_frequencies(timestep,num_steps)
            fs = double(0:(num_steps - 1)) * (pi/(timestep*double(num_steps)));
            fs = fs/16;
        end
        
        function [fs,amps] = randomize_amplitudes(frequencies,spectra,rand_phase)
            amps = zeros(numel(spectra),numel(frequencies));
            fs = reshape(frequencies,1,numel(frequencies));
            
            for i = 1:numel(spectra)
                widths = sqrt(spectra{i}(frequencies));
                norm = sqrt(sum(widths.^2));
                widths = widths * (widths(1) / norm);
                widths(1) = 0;
                amps(i,:) = randn([numel(frequencies),1]) .* reshape(widths,[numel(widths),1]);
            end
            
            %If implementing an amplitude cutoff
            cut = 1.e-9;
            disc = (max(abs(amps),[],1) < cut);
            amps(:,disc) = [];
            fs(disc) = [];
            %******************************
            
            if rand_phase
                phases = exp(1i * rand(size(amps)) * 2 * pi);
                amps = amps .* phases;
            end
            
        end
    end
end

