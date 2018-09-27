classdef TimeEvolution_Noise < TimeEvolution
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        hams
        fourier_amplitudes
        frequencies
        num_channels
        static_ham
        t_phases
        times
    end
    
    methods
        function obj = TimeEvolution_Noise(timestep,num_steps,max_exp,hams,spectra,reals,rand_phase,static_ham)
            obj@TimeEvolution(timestep,num_steps,max_exp);
            assert(iscell(hams)); assert(iscell(spectra));
            assert(numel(hams) == numel(spectra),...
                'Must provide equal numbers of Hamiltonians and spectra');
            obj.num_channels = numel(hams);
            obj.hams = hams;
            obj.static_ham = static_ham;
            obj.frequencies = TimeEvolution_Noise.generate_frequencies(obj.timestep,obj.num_steps);
            obj.fourier_amplitudes = zeros(numel(hams),numel(obj.frequencies),reals);
            for i = 1:numel(spectra)
                amps_i = TimeEvolution_Noise.randomize_amplitudes(...
                    obj.frequencies,spectra{i},reals,rand_phase);
                obj.fourier_amplitudes(i,:,:) = reshape(amps_i,[1,size(amps_i)]);
            end
        end
        
        function ham = calculate_hamiltonian(obj,t_step,real_number)
%             if isempty(real_number)
%                 error('Realization number not given');
%             end
            phases = obj.t_phases(t_step+1,:);
            ham = obj.static_ham;
            for i = 1:obj.num_channels
                signal_i = imag(phases * (obj.fourier_amplitudes(i,:,real_number).'));
                ham = ham + signal_i*obj.hams{i};
            end
        end
        
        function sigs = calculate_signals(obj,t_step,channel_number,real_number)
            phases = obj.t_phases(t_step+1,:);
            sigs = imag(phases * (obj.fourier_amplitudes(channel_number,:,real_number).'));
        end
        
        function obj = allocate_phases(obj,times)
            obj.t_phases = exp(1i * obj.frequencies .* reshape(times,[numel(times),1]));
            obj.times = times;
        end
        
    end
    
    methods (Static)
        function spectra = generate_poissonians(widths,amplitudes)
            assert(numel(widths) == numel(amplitudes),...
                'Must provide equal numbers of widths and amplitudes');
            spectra = cell(size(widths));
            for i = 1:numel(widths)
                spectra{i} = @(x) amplitudes(i)*widths(i)*((1 + x.^2 * widths(i).^2).^(-1));
            end
        end
        
        %Angular frequencies (times 2 pi)
        function fs = generate_frequencies(timestep,num_steps)
            fs = (1:(num_steps - 1)) * (2*pi/(timestep*num_steps));
        end
        
        function amps = randomize_amplitudes(frequencies,spectrum,num_reals,rand_phase)
            widths = sqrt(spectrum(frequencies));
            amps = randn([numel(frequencies),num_reals]) .* ...
                reshape(widths,[numel(widths),1]);
            if rand_phase
                phases = exp(1i * rand(size(amps)) * 2 * pi);
                amps = amps .* phases;
            end
        end
    end
end

