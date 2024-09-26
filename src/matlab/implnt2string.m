function implntstring = implnt2string(implnt)
	switch implnt
		case 0
			implntstring = "ZBC 2009/14 PLA";
		case 1
			implntstring = "True PLA";
		case 2
			implntstring = "Parallel exponential PLA (numerically optimized)";
		case 3
			implntstring = "Parallel exponential PLA (heuristic)";
	end
end

