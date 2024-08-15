%% fig_main_c.m
% Figure depicting performance benefits of PLA approximation schemes.

% Configure stimulus and simulation parameters
durs = 1e-2 .* 2 .^ (0.0:0.5:10.0);   % simulation durations (s)
n_rep = 10;                          % number of repeats
fs = 100e3;                          % sampling rate (Hz)
n_dur = length(durs);
implnts = [0, 1, 2];                 % PLA mode (old, true, new)
n_implnt = length(implnts);

% Pre-allocate storage
runtimes_ihc = zeros(n_dur, n_implnt, n_rep);
runtimes_anf = zeros(n_dur, n_implnt, n_rep);
runtimes_tot = zeros(n_dur, n_implnt, n_rep);

% Loop over durations
for idx_dur = 1:n_dur
	% Loop over models
	for idx_implnt = 1:n_implnt
		% Loop over repeats
		for idx_rep = 1:n_rep
			% Time execution of IHC part
			tic; 
			sim_ihc_zbc2014(zeros(1, round(durs(idx_dur)*fs)), 1000.0); 
			runtimes_ihc(idx_dur, idx_implnt, idx_rep) = toc; 

			% Time execution of ANF part
			tic; 
			sim_an_zbc2023(zeros(1, round(durs(idx_dur)*fs)), 1000.0, implnt=implnts(idx_implnt)); 
			runtimes_anf(idx_dur, idx_implnt, idx_rep) = toc; 

			% Record total execution time
			runtimes_tot(idx_dur, idx_implnt, idx_rep) = ...
				runtimes_ihc(idx_dur, idx_implnt, idx_rep) + ...
				runtimes_anf(idx_dur, idx_implnt, idx_rep);
		end
	end
end

% Create figure
figure;

% Loop through implnts and plot
for idx_implnt = 1:n_implnt
	% Calculate results
	mu = mean(runtimes_anf(:, idx_implnt, :), 3);                 % average runtime across reps
	err = std(runtimes_anf(:, idx_implnt, :), 0, 3)/sqrt(n_rep);  % standard error of average

	% Plot errorbars
	errorbar( ...
		durs, ...
		mu, 1.96 * err, ...
		"Color", implnt2color(implnts(idx_implnt)), ...
		"LineWidth", 3.0 ...
	); hold on;
end

% Plot unity line
plot(10 .^ (-3.0:0.01:3.0), 10 .^ (-3.0:0.01:3.0), 'k'); hold off;  % plot unity line

% Add text labels for comparisons
% mu_true = mean(runtimes_anf(end, idxs_implnts(1), :), 3);
% mu_old = mean(runtimes_anf(end, idxs_implnts(2), :), 3);
% mu_new = mean(runtimes_anf(end, idxs_implnts(3), :), 3);

% Add detailing
grid on;
set(gca, "xscale", "log");
set(gca, "yscale", "log");
ylim([1e-3, 1e2]);
xlim([5e-2, 2e1]);
xlabel("Stimulus duration (s)");
ylabel("Compute time (s)");
xticks([1e-2 1e-1 1e0 1e1]);
yticks([1e-2 1e-1 1e0 1e1]);

% Save to disk
export_at_size(gcf, "figs\fig2c.png", [3.0 2.7], 600);