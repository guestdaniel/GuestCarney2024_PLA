function y = ffGn_rochester(N,tdres,Hinput,noiseType,spont,model_version)
% ffGn_rochester: Fast (exact) fractional Gaussian noise and Brownian motion generator.
%	Y = ffGn_rochester(N, tdres, Hinput, noiseType, spont, model_version)
%   returns a vector containing a sequence of fractional Gaussian
%	noise or fractional Brownian motion generated at 10 samples/sec and
%	upsampled to a sample rate of 1/tdres.  The generation process uses an FFT
%	which makes it very fast.  The input arguments are:
%
%		N			  is the length of the output sequence.
%       tdres         is the time resolution (1/sampling rate)
%		Hinput		  is the "Hurst" index of the resultant noise (0 < H <= 2).
%                     For 0 < H <= 1, the output will be fractional Gaussian
%                     noise with Hurst index H.  For 1 < H <= 2, the output will
%                     be fractional Brownian motion with Hurst index H-1.
%                     Either way, the power spectral density of the output will
%                     be nominally proportional to 1/f^(2H-1).
%		noiseType	  is 0 for fixed fGn noise and 1 for variable fGn 
%		spont		  is the spontaneous rate
%		model_version is the version of sigma to apply, 2014 or 2018
%
% 	References: Davies & Harte (1987); Beran (1994); Bardet et al., 2002
%	This method is based on an embedding of the covariance matrix in a circulant
%	matrix.
%
%   Copyright 2003-2005 by B. Scott Jackson
%   History:
%   Revision: 2.1    29 September 2023 by D. Schwarz and D. Guest.
%                    Changed function name to one that will be widely useful for
%                    all our modeling needs.  Changed variable 'version' to
%                    'model_version' and require it to be either 2014 or 2018.
%   Revision: 2.0    Date: 19 August 2021 by D. Schwarz.
%                    Changed function name.  Changed fifth argument name to
%                    spont (spontaneous rate) and added sixth argument, version,
%                    to switch between two different uses of this function: in
%                    model_Synapse (version 1) and model_Synapse_BEW2018
%                    (version 2).  Moved the persistent statement to its
%                    preferred location.  Now create two different random number
%                    generator (RNG) streams, one with a shuffled seed and the
%                    other with a fixed seed, so that the global RNG is not
%                    affected by using the old V4 algorithm (which is required
%                    here so this function does the same thing today that it did
%                    in 2012).  Used the assert function for argument checking.
%   Revision: 1.4    Date: November 27, 2012 by M. S. A. Zilany : noiseType
%                       has been added
%   Revision: 1.3    Date: Aug 28, 2008 by M. S. A. Zilany
%                    Sigma is defined for diff. sponts (mu) and Resampling has
%                    been introduced to be compatible with the AN model
%   Revision: 1.2    Date: March 14, 2005
%   Rev. 1.2 - 3/14/05 - Added some documentation and input argument checking.
%   Rev. 1.1 - 9/15/04 - Added persistent variables and associated "if" statement.
%   Rev. 1.0 - 2/11/03 - Original version.


persistent Zmag Nfft Nlast Hlast rs_shuffled rs_fixed

%---- Check input arguments ---------- %

narginchk(6,6)

if ~isscalar(N) || ~isscalar(Hinput) || ~isnumeric(N) || ~isnumeric(Hinput) ...
		|| ~isreal(N) || ~isreal(Hinput) || ~isfinite(N) || ~isfinite(Hinput)
	error('All input arguments must be finite real scalars.')
end

assert(N > 0,'Length of the return vector must be positive.')
assert(tdres <= 1,'Original sampling rate should be checked.')
assert(Hinput > 0 && Hinput <= 2,'The Hurst parameter must be in the interval (0,2].')
assert(ismember(noiseType,0:1),'The noiseType must be 0 or 1.')
assert(spont > 0,'The spontaneous rate must be positive.')
assert(ismember(model_version,[2014 2018]),'Model version must be 2014 or 2018.')


% Downsampling number of points to match those of Scott Jackson (Tj = 0.1s).
Tj = 0.1; % s
resamp = ceil(Tj/tdres);
nop = N;
N = max(ceil(N/resamp) + 1,10);

% Determine whether fGn or fBn should be produced.
if Hinput <= 1
	H = Hinput;
	is_fBn = false;
else
	H = Hinput - 1;
	is_fBn = true;
end

% Create random number streams the first time this function is called.
if isempty(rs_shuffled)
	rs_shuffled = RandStream('mt19937ar','Seed','shuffle');
	rs_fixed = RandStream('v4');
end

% Select random number generator.
if noiseType == 0
	% Use V4 random number generator with a fixed seed (37).
	rs_fixed.reset(37)
	rs = rs_fixed;
elseif noiseType == 1
	% Use mt19937ar (default) random number generator with a shuffled seed.
	rs = rs_shuffled;
end

% Calculate the fGn.
if H == 0.5
	% If H = 0.5, then fGn is equivalent to white Gaussian noise.
	y = rs.randn(1,N);
else
	% If this function was already in memory before being called this time, AND
	% the values for N and H are the same as the last time it was called, then
	% the following (persistent) variables do not need to be recalculated.  This
	% was done to improve the speed of this function, especially when many
	% samples of a single fGn (or fBn) process are needed by the calling
	% function.
	if isempty(Hlast) || N ~= Nlast || H ~= Hlast
		% The persistent variables must be (re-)calculated.
		Nfft = 2.^nextpow2(2*(N-1));
		NfftHalf = Nfft/2;
		
		k = [0:NfftHalf, (NfftHalf-1):-1:1];
		Zmag = 0.5 * ( (k+1).^(2.*H) - 2.*k.^(2.*H) + (abs(k-1)).^(2.*H) );
		Zmag = real(fft(Zmag));
		if any(Zmag < 0)
			error('The FFT of the circulant covariance had negative values.');
		end
		Zmag = sqrt(Zmag);
		
		% Store N and H values in persistent variables for use during subsequent
		% calls to this function.
		Nlast = N;
		Hlast = H;
	end
	
	Z = Zmag.*complex(rs.randn(1,Nfft),rs.randn(1,Nfft));
	
	y = sqrt(Nfft)*real(ifft(Z));
	y = y(1:N);
end

% Convert the fGn to fBn, if necessary.
if is_fBn
	y = cumsum(y);
end

% Resampling back to original (1/tdres) to match with the AN model
y = resample(y,resamp,1);

% Set the standard deviation, sigma.
if model_version == 2014 % 2014 model
	if spont < 0.5
		sigma = 3; % 5
	elseif spont < 18
		sigma = 30; % 50   % 7 when added after powerlaw
	else
		sigma = 200; % 40 when added after powerlaw
	end
elseif model_version == 2018 % BEZ2018 model
	if spont < 0.2
		sigma = 1; % 5
	elseif spont < 20
		sigma = 10;
	else
		sigma = spont/2;
	end
else
	error('Unknown model version.')
end

% Apply sigma and clip unneeded samples.
y = sigma*y(1:nop);
