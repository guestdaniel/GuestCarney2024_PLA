function stim = cbc_roving_samt(args)
% CBC_ROVING_SAMT(args) Synthesizes a novel stimulus used to probe
% power-law adaptation in the Zilany, Bruce, and Carney (2014) model and
% the Guest and Carney (2023) power-law adaptation paper. 
%
% The stimulus is a variation on an amplitude-modulated tone whereby the
% the sound level of waveform corresponding to each modulation cycle is
% randomized over a uniform range and the time between two successive
% modulation cycles is also randomized over a uniform range instead of
% being fixed at 0 s (as in traditional SAM).
	arguments
		args.f = 1000.0;              % carrier frequency (Hz)
		args.fm = 10.0;               % modulation frequency (Hz)
		args.m = 1.0;                 % modulation depth
		args.dur = 1.0;               % duration (s)
		args.fs = 100e3;              % sampling rate (Hz)
		args.rove_range = [20 70]     % min/max rove (dB SPL)
		args.isi_range = [0.0 0.1];   % isi range (s)
	end

	% Create time vectors for single cycle of modulator and entire tone
	t_one_cycle = 0.0:(1/args.fs):((1/args.fm) - 1/args.fs);
	t = 0.0:(1/args.fs):(args.dur - 1/args.fs);

	% Synthesize the carrier
	carrier = sin(2*pi * args.f * t);

	% Create the stimulus in steps of one mod cycle + one ISI
	stim = zeros(size(carrier));
	t_curr = 0.0;
	while t_curr < args.dur
		% Select ISI randomly
		isi = (args.isi_range(2) - args.isi_range(1)) * rand + args.isi_range(1);

		% Select the sound level
		level = (args.rove_range(2) - args.rove_range(1)) * rand + args.rove_range(1);

		% Create single cycle of modulator and ISI
		modulator_one_cycle = [sin(2*pi * t_one_cycle * args.fm - pi/2) -1*ones(round(isi*args.fs), 1)'];

		% Calculate indices 
		idx_start = round(t_curr*args.fs) + 1;
		idx_end = idx_start + length(modulator_one_cycle) - 1;
		idxs_valid = idx_start:min(length(t), idx_end);

		% Create stimulus (unscaled) for single mod cycle + isi
		stim_this_cycle = (1 + args.m*modulator_one_cycle(1:length(idxs_valid))) .* carrier(idxs_valid);

		% Scale and store stimulus in output vector
		stim_this_cycle = 20e-6 * 10^(level/20.0) * stim_this_cycle/rms(stim_this_cycle);
		stim(idxs_valid) = stim_this_cycle;

		% Increment time
		t_curr = idx_end/args.fs;
	end
end