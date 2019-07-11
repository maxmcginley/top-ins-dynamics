classdef TimeEvolution_Noise < TimeEvolution
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hams
        fourier_amplitudes
        frequencies
        num_channels
        static_ham
        signal
        times
    end
    
    methods
        function obj = TimeEvolution_Noise(timestep,num_steps,max_exp,hams,spectra,rand_phase,static_ham)
            obj@TimeEvolution(timestep,num_steps,max_exp);
            assert(iscell(hams)); assert(iscell(spectra));
            assert(numel(hams) == numel(spectra),...
                'Must provide equal numbers of Hamiltonians and spectra');
            obj.num_channels = numel(hams);
            obj.hams = hams;
            obj.static_ham = static_ham;
            fqs = TimeEvolution_Noise.generate_frequencies(obj.timestep,obj.num_steps);
            [obj.frequencies, obj.fourier_amplitudes] = TimeEvolution_Noise.randomize_amplitudes(...
                    fqs,spectra,rand_phase);
        end
        
        function ham = calculate_hamiltonian(obj,t_step,~)
%             if isempty(real_number)
%                 error('Realization number not given');
%             end
            ham = obj.static_ham;
            
            for i = 1:obj.num_channels
                %signal_i = imag(obj.signal(i,mod(t_step,size(obj.signal,2))+1));
                signal_i = sum(imag(exp(1i * t_step * obj.timestep * obj.frequencies) .* obj.fourier_amplitudes));
                ham = ham + signal_i*obj.hams{i};
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
        
        function amp = spec_function(x,wid)
            %amp = (2*wid/pi)*((wid.^2 + x.^2).^(-1));
            amp = (1/(2*wid)) * exp(-2*(wid^(-1))*abs(x));
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
        end
        
        function [fs,amps] = randomize_amplitudes(frequencies,spectra,rand_phase)
            amps = zeros(numel(spectra),numel(frequencies));
            fs = reshape(frequencies,1,numel(frequencies));
            
            for i = 1:numel(spectra)
                widths = sqrt(spectra{i}(frequencies));
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

