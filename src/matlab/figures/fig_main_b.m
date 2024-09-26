%% fig_main_b.m
% Figure depicting how a running-error metric compares for old versus new
% power-law adaptation approximation

% Configure stimulus and simulation parameters
dur = 100.0;                           % stim duration (s)
dur_post = 0.25;                       % extra sim dur after stim offset (s)
cf = 10e3;                             % tone frequency (Hz)
fibertype = 3;                         % HSR fiber
species = 1;                           % cat model
noisetype = 0;		                   % no fGn
implnts = [0, 1, 2, 3];                % PLA mode (old, orig, new-opt, new-heur)
n_implnt = length(implnts);
n_rep = 1;                             % how many different waveforms to simulate
fs = 100e3;

% Configure analysis windowing parameters
dur_window = 0.2;  % seconds
len_window = round(dur_window*100e3);
overlap = 0.5;     % proportion
dur_skip = dur_window * overlap;
len_skip = round(dur_skip*100e3);

% Construct time axis (will need later)
t = 0.0:(1/100e3):((dur+dur_post)-1/100e3);

% Calculate window times and indices
idxs = 1:len_skip:(length(t)-len_window);
t_windowed = ((1:len_skip:(length(t)-len_window)) - 1) / 100e3;

% Pre-allocate result storage
errors_windowed = zeros(length(t_windowed), 3, n_rep);

% Do main loop
for idx_rep = 1:n_rep
	% Pre-allocate storage for response waveforms
	responses = cell(n_implnt, 1);
	
	% Synthesize CBC roving SAMT
	%stim = cbc_roving_samt(f=cf, dur=dur, fm=10.0, m=1.0);
	stim = sam_tone(cf, 10.0, 1.0, 50.0, dur, 0.01, fs);
	
	% Simulate responses at IHC stage
	ihc_orig = sim_ihc_zbc2014(stim, cf, species=species, dur=(dur+dur_post));
	
	% Loop through implnts and sim each AN response
	for idx_implnt = 1:n_implnt
		responses{idx_implnt} = sim_an_zbc2023( ...
			ihc_orig, ...
			cf, ...
			... Parameters below should yield 2014 model w/o fGn, but
			... different PLA depending on parameter value
			fibertype=fibertype, ...
			noisetype=noisetype, ...
			implnt=implnts(idx_implnt) ...
		);
	end
	
	% Compute error for each implementation
	errors = zeros(length(t), 3);
	errors(:, 1) = responses{2} - responses{1};  % error for old approximate PLA
	errors(:, 2) = responses{2} - responses{3};  % error for new approximate PLA (optimized)
	errors(:, 3) = responses{2} - responses{4};  % error for new approximate PLA (heuristic)

	% Compute windowed MAE
	idx_windowed = 1;
	temp = zeros(length(t_windowed), 3);
	for idx = 1:len_skip:(length(t)-len_window)
		temp(idx_windowed, 1) = rms(errors(idx:(idx+len_window), 1));
		temp(idx_windowed, 2) = rms(errors(idx:(idx+len_window), 2));
		temp(idx_windowed, 3) = rms(errors(idx:(idx+len_window), 3));
		idx_windowed = idx_windowed + 1;
	end

	% Store results
	errors_windowed(:, :, idx_rep) = temp;
end

% Create figure
figure;

% % Plot each individual trace as a thin partially transparent line
% for idx_rep = 1:n_rep
% 	plot(t_windowed, errors_windowed(:, 1, idx_rep), color=[hex2rgb(char(implnt2color(0))) 0.05], linewidth=1.0); hold on;
% 	plot(t_windowed, errors_windowed(:, 2, idx_rep), color=[hex2rgb(char(implnt2color(2))) 0.05], linewidth=1.0);
% 	plot(t_windowed, errors_windowed(:, 3, idx_rep), color=[hex2rgb(char(implnt2color(2))) 0.05], linewidth=1.0);
% end

% Plot mean as thicker solid line
plot(t_windowed, mean(errors_windowed(:, 1, :), 3), color=implnt2color(0), linewidth=2.0); hold on;
plot(t_windowed, mean(errors_windowed(:, 2, :), 3), color=implnt2color(2), linewidth=2.0);
plot(t_windowed, mean(errors_windowed(:, 3, :), 3), color=implnt2color(3), linewidth=2.0, linestyle=implnt2linestyle(3));

% Adjust plot to log time axis
set(gca, "Xscale", "log");

% Add labels
grid on;
xlabel("Time (s)");
ylabel("RMS error (sp/s)");
xlim([0 dur]);

% Save to disk
export_at_size(gcf, "figs\fig2b.png", [3, 2.7], 600);