function [rate, var, spikes] = sim_an_zbc2023(x, cf, args)
% SIM_AN_ZBC2023(...) Simulates ANF rates/spikes using auditory-periphery 
% model of Zilany, Bruce, and Carney (2014), optionally using the 
% Guest and Carney (2023) power-law adaptation approximation. Returns
% outputs as column vectors.
%
% [rate, var, spikes] = sim_an_zbc2023(...) returns the outputs of the ANF
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
		args.fs_synapse (1,1) double = 10e3
		args.tau_slow = 5e-4 .* 10 .^ (-1:(1/exp(1)):5);
		args.w_slow = [220.75631209760925 155.36035581822114 725.4659020418738 286.7067395494722 554.5094952969946 200.1877115289392 77.38005673143523 64.40005070188093 7.367301456151062 12.804999507365363 2.559156443618358 1.4833275422533558 0.6629878278644696 0.3614808470276771 0.08264813156030201 0.003229786278458085 0.07977841451420807]
   		args.tau_fast = 1e-1 .* 10 .^ (-1:(1/exp(1)):5);
		args.w_fast = [220.75631209760925 155.36035581822114 725.4659020418738 286.7067395494722 554.5094952969946 200.1877115289392 77.38005673143523 64.40005070188093 7.367301456151062 12.804999507365363 2.559156443618358 1.4833275422533558 0.6629878278644696 0.3614808470276771 0.08264813156030201 0.003229786278458085 0.07977841451420807]
	end

	% Pass inputs to the Mex wrapper, model_Syanapse_2023
    [rate, var, spikes] = model_Synapse_2023(...
		x', ...
		cf, ...
		args.nrep, ...
		1/args.fs, ...
		args.fibertype, ...
		args.noisetype, ...
		args.implnt, ...
		args.fs_synapse, ...
		args.tau_slow, ...
		args.w_slow, ...
		args.tau_fast, ...
		args.w_fast, ...
		length(args.tau_slow) ...
	);

	% Transform row-vector outputs into column-vector outputs
	rate = rate';
	var = var';
	spikes = spikes';
end

