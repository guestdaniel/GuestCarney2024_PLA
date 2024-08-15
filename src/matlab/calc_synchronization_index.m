function R = calc_synchronization_index(spk_times, f)
	R = abs(1/length(spk_times) * sum(exp(1i * 2*pi * f .* spk_times)));
end

