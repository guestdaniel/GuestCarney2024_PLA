%% fig_main_a.m
% Figure depicting how new PLA adaptation scheme better handles responses
% to a challenging stimulus at both short and long time scales.

function fig = fig_main_a()

% Configure stimulus and simulation parameters
dur = 2.5;                             % stim duration (s)
dur_post = 0.25;                       % extra sim dur after stim offset (s)
cf = 10e3;                             % tone frequency (Hz)
fibertype = 3;                         % HSR fiber
species = 1;                           % cat model
noisetype = 0;		                   % no fGn
implnts = [0, 1, 2];                   % PLA mode (old, true, new)
n_implnt = length(implnts);
fs = 100e3;                            % sampling rate (Hz)

% Pre-allocate storage for response waveforms
responses = cell(n_implnt);

% Synthesize CBC roving SAMT
stim = cbc_roving_samt(f=cf, dur=dur, fm=10.0, m=1.0);

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

% Create time axis
t = 0.0:(1/100e3):((dur+dur_post)-1/100e3);

% Create and configure figure
dur_skip = 0.01;  % how much of a gap to leave between time skip
fig = figure;
tiledlayout(5, 1, "TileSpacing", "none", "Padding", "tight");

% Set up axes
ax_stim = nexttile;
hold(ax_stim, "on");
set(ax_stim, "Visible", "off");
ax_sim = nexttile(2, [4, 1]);
hold(ax_sim, "on");

% We want to plot the initial ~1/4 of the response and the final ~1/4 of
% the response on a single time axis, so we'll place a time-skip in the
% axis around the 600ms or so. Here, we determine indices for first 1/4 
% and last 1/4
idxs_init = t < (t(end)/4);
idxs_final = t > (t(end)*3/4);

% Plot initial part (stimulus)
plot( ...
	ax_stim, ...
	t(idxs_init), ...
	stim(idxs_init), ...
	color="black" ...
);

% Plot initial part (stimulus)
t_final = t(idxs_init) + dur_skip + t(end)/4;
plot( ...
	ax_stim, ...
	t_final(1:length(stim(idxs_final(1:length(stim))))), ...
	stim(idxs_final(1:length(stim))), ...
	color="black" ...
);

% Loop through implnts and plot each result
for idx_implnt = 1:n_implnt
	% Plot initial part (simulations)
	plot( ...
		ax_sim, ...
		t(idxs_init), ...
		responses{idx_implnt}(idxs_init), ...
		linestyle=implnt2linestyle(implnts(idx_implnt)), ...
		linewidth=implnt2linewidth(implnts(idx_implnt)), ...
		color=implnt2color(implnts(idx_implnt)) ...
	);

	% Plot final part on "cheater's" time axis with time skip (simulations)
	plot( ...
		ax_sim, ...
		t(idxs_init) + dur_skip + t(end)/4, ...
		responses{idx_implnt}(idxs_final), ...
		"HandleVisibility", "off", ...
		linestyle=implnt2linestyle(implnts(idx_implnt)), ...
		linewidth=implnt2linewidth(implnts(idx_implnt)), ...
		color=implnt2color(implnts(idx_implnt)) ...
	);
end

% Manually adjust xticks to get real and skipped time ticks
xticks_init = 0.0:0.1:(t(end)/4);
xticks_final = 2.1:0.1:t(end);
xticks_final_adj = xticks_final + dur_skip - t(end)*2/4;
set(ax_stim, "Xtick", [xticks_init xticks_final_adj]);
set(ax_stim, "Xticklabel", [string(xticks_init) string(xticks_final)])
set(ax_sim, "Xtick", [xticks_init xticks_final_adj]);
set(ax_sim, "Xticklabel", [string(xticks_init) string(xticks_final)])

% Set labels and such
xlabel(ax_sim, "Time (s)");
ylabel(ax_sim, "Firing rate (sp/s)");

% Set limits
set(ax_stim, "Xlim", [0 1.35]);
set(ax_sim, "Xlim", [0 1.35]);

% Add legend
legend(["Old approximate PLA", "True PLA", "Parallel exponential PLA"]);

export_at_size(gcf, "figs\fig2a.png", [7, 4], 600);
end