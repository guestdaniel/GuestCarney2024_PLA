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
	end

	% Determine tau and w for PLA approximation based on requested implnt
	if args.implnt == 0 
		% If we pass 0 (old approximation) or 1 (true PLA), simply set to 0
		tau_slow = [0.0]; w_slow = [0.0];
		tau_fast = [0.0]; w_fast = [0.0];
		implnt = 0;
	elseif args.implnt == 1
		% If we pass 0 (old approximation) or 1 (true PLA), simply set to 0
		tau_slow = [0.0]; w_slow = [0.0];
		tau_fast = [0.0]; w_fast = [0.0];
		implnt = 1;
	elseif args.implnt == 2
% 		tau_slow = 5e-4 .* 10 .^ (-1:(1/exp(1)):5);
% 		w_slow = [220.75631209760925 155.36035581822114 725.4659020418738 286.7067395494722 554.5094952969946 200.1877115289392 77.38005673143523 64.40005070188093 7.367301456151062 12.804999507365363 2.559156443618358 1.4833275422533558 0.6629878278644696 0.3614808470276771 0.08264813156030201 0.003229786278458085 0.07977841451420807];
%    		tau_fast = 1e-1 .* 10 .^ (-1:(1/exp(1)):5);
% 		w_fast = [220.75631209760925 155.36035581822114 725.4659020418738 286.7067395494722 554.5094952969946 200.1877115289392 77.38005673143523 64.40005070188093 7.367301456151062 12.804999507365363 2.559156443618358 1.4833275422533558 0.6629878278644696 0.3614808470276771 0.08264813156030201 0.003229786278458085 0.07977841451420807];
% 		implnt = 2;
 		tau_slow = 5e-4 .* 10 .^ (0:(1/exp(1)):5);
		w_slow = [1054.1349144510866 235.42021095822022 351.3091743124357 99.00123234954474 55.18423650003196 28.99454378212968 6.556134147763605 6.558380224204848 1.1576087874250394 0.995488845827021 0.3588672871386332 0.1573449044190812 0.010428823220777147 0.08773889583510958];
		tau_fast = 1e-1 .* 10 .^ (0:(1/exp(1)):5);
		w_fast = [6.106637716398411 1.1558964083697898 1.3095958543425545 0.785695677692722 0.21835528692662018 0.10344785429373701 0.08413927488982781 0.001596356536824024 0.018886711336962816 0.0008089617759213521 0.002806098243601203 0.0006529927704604911 3.13953727422695e-5 0.0004490670084957763];
		implnt = 2;
	elseif args.implnt == 3
		tau_slow = 5e-4 .* 10 .^ (0:(1/exp(1)):5);
		w_slow = 1 ./ (tau_slow + 5e-4);
		c =  2*5e-4 * sum(w_slow .* exp(-5e-4 ./ tau_slow));
		w_slow = w_slow ./ c;
		tau_fast = 1e-1 .* 10 .^ (0:(1/exp(1)):5);
		w_fast = 1 ./ (tau_fast + 1e-1);
		c =  2*1e-1 * sum(w_fast .* exp(1e-1 ./ tau_fast));
		w_fast = w_fast ./ c;
		implnt = 2;
	end

	% Pass inputs to the Mex wrapper, model_Syanapse_2023
    [rate, var, spikes] = model_Synapse_2023(...
		x', ...
		cf, ...
		args.nrep, ...
		1/args.fs, ...
		args.fibertype, ...
		args.noisetype, ...
		implnt, ...
		args.fs_synapse, ...
		tau_slow, ...
		w_slow, ...
		tau_fast, ...
		w_fast, ...
		length(tau_slow) ...
	);

	% Transform row-vector outputs into column-vector outputs
	rate = rate';
	var = var';
	spikes = spikes';
end

