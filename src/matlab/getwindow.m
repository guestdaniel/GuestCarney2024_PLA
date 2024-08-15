function [t_low, t_high] = getwindow(cf, dur, level, type)
% GETWINDOW Get time points for response windows as function of cf for 
% stimulus of duration dur and overall level level. Suitable only for CFs 
% below 16 kHz and levels below 80 dB SPL.
%
%	[t_low, t_high] = GETWINDOW(cf, dur, "onset") attempts to return the
%	time range corresponding to the onset response at given CF 
%
%	[t_low, t_high] = GETWINDOW(cf, dur, "ss") attempts to return the time
%	range corresponding to the steady-state response at given CF
%
% Steady-state windows are defined as everything from 10-ms and beyond
%
% Onset windows are more complex. First, we define a t_low as a function of
% sound level and CF. Specifically, t_low is set to the first-spike
% latency (as a function of CF) plus a linear function of sound level that
% further delays t_low for at low sound levels. t_high is set at t_low plus
% a 2 ms window afterward plus an additional 0.8 ms window per octaves below
% 16 kHz.

	% First, determine expected first-spike latency based on cat AN data
	A0 = 3.0;
	A1 = 12.5;
	x = 11.9 * log10(0.8 + cf / 456);
	delay = A0 * exp(-x/A1) * 1e-3;

	% Next, based on type of window, determine start and stop times
	switch type
		case "onset"
			% Set lower side of window to delay time
			t_low = delay + max(0, 80-level) * 1e-5;

			% Set up upper side as function of CF
			t_high = t_low + 2e-3 + max(0, log2(16e3/cf)) * 0.8e-3;
		case "ss"
			t_low = 0.01;
			t_high = dur;
		case "fullresponse"
			t_low = 0.0;
			t_high = dur + 0.1;
	end
	
end

