function [rate, var, spikes] = sim_an_zbc2025(x, cf, args)
% SIM_AN_ZBC2025(...) Simulates ANF rates/spikes using auditory-periphery 
% model of Zilany, Bruce, and Carney (2014), optionally using the 
% Guest and Carney (2024) power-law adaptation approximation. Returns
% outputs as column vectors.
%
% [rate, var, spikes] = sim_an_zbc2025(...) returns the outputs of the ANF
% stage of the model, which are (respectively) the waveform of
% instantaneous spike rates, the waveform of instaneous spike variances,
% and a simulated peristimulus time histogram (PSTH). 
%
% Arguments:
% - x: inner hair cell potential (a.u.), size (1, totalstim)
% - cf: characteristic frequency (Hz)
% - args.nrep: how many times to run the simulation
% - args.fs: sampling rate (Hz)
% - args.fiber_type: spontaneous-rate type, either 1 (LSR), 2 (MSR), or 
%		3 (HSR)
% - args.noisetype: what type of noise to use in the synapse stage, either 
%		frozen fractional Gaussian noise (0) or fresh fractional Gaussian 
%		noise (1)
% - implnt: how to implement power-law adaptation, either original
%		approximate (0), true power-law adaptation (1), or new 
%		approximate (2). New approximate is the recommended setting for 
%		most applications.
    arguments
        x (:, 1) double 
        cf (1,1) double
        args.nrep (1,1) double = 1
        args.fs (1,1) double = 100e3
        args.fibertype (1,1) double = 3 
        args.noisetype (1,1) double = 0 
        args.implnt (1,1) double = 2
	end

	% Pass inputs to the Mex wrapper, model_Syanapse_2023
    [rate, var, spikes] = model_Synapse_v2025a(...
		x', ...
		cf, ...
		args.nrep, ...
		1/args.fs, ...
		args.fibertype, ...
		args.noisetype, ...
		args.implnt ...
	);

	% Transform row-vector outputs into column-vector outputs
	rate = rate';
	var = var';
	spikes = spikes';
end

