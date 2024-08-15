function signal_out = raised_cosine_ramp(signal, dur_ramp, fs)
% RAISED_COSINE_RAMP Ramps an input signal with raised-cosine onset/offset ramps
%
% signal_out = RAISED_COSINE_RAMP(signal, dur_ramp, fs) adds raised-cosine onset/
% offset ramps of duration dur_ramp to a signal with sampling rate fs.
%
    t = linspace(0.0, 0.25-1/fs, round(fs*dur_ramp));
    ramp_segment = sin(2*pi*t).^2;
    ramp = [ramp_segment ones(1, length(signal) - length(ramp_segment)*2) fliplr(ramp_segment)];
    signal_out = ramp .* signal;
end

