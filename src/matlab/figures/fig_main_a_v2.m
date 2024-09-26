%% fig_main_a_v2.m
% Figure depicting how new PLA adaptation scheme better handles responses
% to a challenging stimulus at both short and long time scales.

% Configure stimulus and simulation parameters
dur = 10.0;                            % stim duration (s)
dur_post = 0.25;                       % extra sim dur after stim offset (s)
cf = 10e3;                             % tone frequency (Hz)
fibertype = 3;                         % HSR fiber
species = 1;                           % cat model
noisetype = 0;		                   % no fGn
implnts = [0, 1, 2, 3];                % PLA mode (old, true, new)
n_implnt = length(implnts);
fs = 100e3;

% Synthesize CBC roving SAMT
stimuli = {...
	sam_tone(cf, 64.0, 1.0, 50.0, dur, 0.01, fs), ...
	scale_dbspl(cosine_ramp(pure_tone(cf, 0.0, dur, fs), 0.01, fs), 50.0), ...
	scale_dbspl(cosine_ramp(randn(dur*fs, 1), 0.01, fs), 50.0) ...
};

% Pre-allocate storage for response waveforms
responses = cell(length(stimuli), n_implnt);

% Simulate responses at IHC stage
ihc_orig = cellfun( ...
	@(x) sim_ihc_zbc2014(x, cf, species=species, dur=(dur+dur_post)), ...
	stimuli, ...
	'UniformOutput', false ...
);

% Loop through implnts and sim each AN response
for idx_stim = 1:length(stimuli)
	for idx_implnt = 1:n_implnt
		responses{idx_stim, idx_implnt} = sim_an_zbc2023( ...
			ihc_orig{idx_stim}, ...
			cf, ...
			... Parameters below should yield 2014 model w/o fGn, but
			... different PLA depending on parameter value
			fibertype=fibertype, ...
			noisetype=noisetype, ...
			implnt=implnts(idx_implnt) ...
		);
	end
end

% Create time axis
t = 0.0:(1/100e3):((dur+dur_post)-1/100e3);

% Create and configure figure
dur_win = 0.05;
wins = {[0.0 0.0 + dur_win], [9.0, 9.0 + dur_win]};
dur_skip = 0.01;  % how much of a gap to leave between time skip
fig = figure;
tl = tiledlayout(5*length(stimuli), 1, "TileSpacing", "none", "Padding", "tight");

% Loop over stimuli and plot each waveform and response
for idx_stim = 1:length(stimuli)
	% Set up axes
	ax_stim = nexttile;
	hold(ax_stim, "on");
	set(ax_stim, "Visible", "off");
	ax_sim = nexttile([3, 1]);
	hold(ax_sim, "on");
	ax_empty = nexttile;
	set(ax_empty, "Visible", "off");

	% Create time axis
	t = 0.0:(1/fs):(length(stimuli{idx_stim})/fs - 1/fs);

	% Loop over windows
	for idx_win = 1:length(wins)
		% Determine indices contained within windows
		idxs = (t >= wins{idx_win}(1)) & (t < wins{idx_win}(2));

		% Create fake time vector for axis
		timevec = t(idxs) + (dur_win + dur_skip)*(idx_win-1) - wins{idx_win}(1);

		% Plot stimuli on top and response on bottom
		plot(ax_stim, timevec, stimuli{idx_stim}(idxs), color="k"); 
		for idx_implnt=1:length(implnts)
			plot( ...
				ax_sim, ...
				timevec, ...
				responses{idx_stim, idx_implnt}(idxs), ...
				linestyle=implnt2linestyle(implnts(idx_implnt)), ...
				linewidth=implnt2linewidth(implnts(idx_implnt)), ...
				color=implnt2color(implnts(idx_implnt)) ...
			)
		end

		% Manually place xticks
		ticks_1 = linspace(0.0, dur_win, 3);
		ticks_2 = ticks_1 + dur_skip + dur_win;
		xticks(ax_sim, [ticks_1 ticks_2])
		xticklabels(ax_sim, [ticks_1 (ticks_2 + wins{2}(1) - dur_win - dur_skip)])
		if idx_stim < length(stimuli)
			xticklabels(ax_sim, []);
		end
	end
end

% Set labels and such
% xlabel(tl, "Time (s)");
% ylabel(tl, "Firing rate (sp/s)");

export_at_size(gcf, "figs\fig3d.png", [4.5, 2.5], 600);