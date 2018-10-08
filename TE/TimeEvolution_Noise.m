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
            obj.frequencies = TimeEvolution_Noise.generate_frequencies(obj.timestep,obj.num_steps);
            obj.fourier_amplitudes = zeros(numel(hams),numel(obj.frequencies));
            for i = 1:numel(spectra)
                amps_i = TimeEvolution_Noise.randomize_amplitudes(...
                    obj.frequencies,spectra{i},rand_phase);
                obj.fourier_amplitudes(i,:) = reshape(amps_i,[1,numel(amps_i)]);
            end
        end
        
        function ham = calculate_hamiltonian(obj,t_step,~)
%             if isempty(real_number)
%                 error('Realization number not given');
%             end
            ham = obj.static_ham;
            for i = 1:obj.num_channels
                signal_i = imag(obj.signal(i,mod(t_step,size(obj.signal,2))+1));
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
        function spectra = generate_poissonians(widths,amplitudes)
            assert(numel(widths) == numel(amplitudes),...
                'Must provide equal numbers of widths and amplitudes');
            spectra = cell(size(widths));
            for i = 1:numel(widths)
                spectra{i} = @(x) (amplitudes(i).^2)*(2*widths(i)/pi)*((widths(i).^2 + x.^2).^(-1));
            end
        end
        
        %Angular frequencies (times 2 pi)
        function fs = generate_frequencies(timestep,num_steps)
            fs = double(0:(num_steps - 1)) * (pi/(timestep*double(num_steps)));
        end
        
        function amps = randomize_amplitudes(frequencies,spectrum,rand_phase)
            widths = sqrt(spectrum(frequencies));
            widths(1) = 0;
            amps = randn([numel(frequencies),1]) .* ...
                reshape(widths,[numel(widths),1]) * (sqrt(2*max(frequencies)/numel(frequencies)));
            if rand_phase
                phases = exp(1i * rand(size(amps)) * 2 * pi);
                amps = amps .* phases;
            end
        end
    end
end

