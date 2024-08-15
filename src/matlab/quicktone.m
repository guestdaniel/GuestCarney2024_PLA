function stim = quicktone(freq, dur, dur_ramp, level, fs)
	arguments
		freq
		dur=0.1
		dur_ramp=0.01
		level=25.0
		fs=100e3
	end
	% Synthesize pure tone with raised-cosine ramp and scale
	t = 0.0:(1/fs):(dur - 1/fs);                 % sample times (s) 
	stim = sin(2*pi * freq * t);
	stim = raised_cosine_ramp(stim, dur_ramp, fs);
	stim = 20e-6 * 10^(level/20.0) * stim/rms(stim);
end

