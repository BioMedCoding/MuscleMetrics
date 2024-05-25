function c = xcorr(a, b, option)
%XCORR	Cross-correlation function estimates.
%	XCORR(A,B), where A and B are length M vectors, returns the
%	length 2*M-1 cross-correlation sequence in a column vector.
%	XCORR(A), when A is a vector, is the auto-correlation sequence.
%	XCORR(A), when A is an M-by-N matrix, is a large matrix with
%	2*M-1 rows whose N^2 columns contain the cross-correlation
%	sequences for all combinations of the columns of A.
%	The zeroth lag of the output correlation is in the middle of the 
%	sequence, at element or row M.
%	By default, XCORR computes a raw correlation with no normalization.
%	XCORR(A,'biased') or XCORR(A,B,'biased') returns the "biased"
%	estimate of the cross-correlation function.  The biased estimate
%	scales the raw cross-correlation by 1/M.
%	XCORR(...,'unbiased') returns the "unbiased" estimate of the
%	cross-correlation function.  The unbiased estimate scales the raw
%	correlation by 1/(M-abs(k)), where k is the index into the result.
%	XCORR(...,'coeff') normalizes the sequence so that the
%	correlations at zero lag are identically 1.0.
%	See also XCOV, CORRCOEF, CONV and XCORR2.

%	Author(s): L. Shure, 1-9-88
%		   L. Shure, 4-13-92, revised
%	Copyright (c) 1984-94 by The MathWorks, Inc.
%	$Revision: 1.7 $  $Date: 1994/01/25 18:00:07 $

%	References:
%	  [1] J.S. Bendat and A.G. Piersol, "Random Data:
%	      Analysis and Measurement Procedures", John Wiley
%	      and Sons, 1971, p.332.
%	  [2] A.V. Oppenheim and R.W. Schafer, Digital Signal 
%	      Processing, Prentice-Hall, 1975, pg 539.

onearray = 1;
if nargin == 1
	option = 'none';
	if min(size(a)) == 1	% a is a vector
		a = [a(:) a(:)];
	else
		onearray = 0;
	end
elseif nargin == 2
	if isstr(b)
		option = b; clear b
		na = max(size(a));
		if min(size(a)) == 1	% a is a vector
			a = [a(:) a(:)];
		else	% a is a matrix
			onearray = 0;
			[m,n] = size(a);
		end
    else	% b is truly a second arg
		if min(size(a)) ~= 1 & min(size(b)) ~= 1
			error('You may only specify 2 vector arrays.')
		end
		option = 'none';
		onearray = 2;
	end
else
	if max(size(a)) ~= max(size(b)) & ~strcmp(option,'none')
		error('OPTION must be ''none'' for different length vectors A and B')
	end
	onearray = 2;
end
% check validity of option
nopt = nan;
if strcmp(option, 'none')
	nopt = 0;
elseif strcmp(option, 'coeff')
	nopt = 1;
elseif strcmp(option, 'biased')
	nopt = 2;
elseif strcmp(option, 'unbiased')
	nopt = 3;
end
if isnan(nopt)
	error('Unknown OPTION')
end
if onearray == 2
	[ar,ac] = size(a);
	na = max([ar ac]);
	nb = max(size(b));
	if na > nb
		b(na) = 0;
	elseif na < nb
		a(nb) = 0;
	end
	a = [a(:) b(:)];
end
[nr, nc] = size(a);
nsq  = nc^2;
mr = 2 * nr - 1;
c = zeros(mr,nsq);
ci = zeros(1,nsq);
cj = ci;
nfft = 2^nextpow2(2*nr);
for i = 1:nc
	atmpi = a(:,i);
	if ~any(any(atmpi))
		real1 = 1;
	else
		real1 = 0;
	end
	atmpi = fft([atmpi(:); zeros(nfft-nr,1)]);
	for j = 1:i
		col = (i-1)*nc+j;
		colaux = (j-1)*nc+i;
		tmp = fft([a(:,j); zeros(nfft-nr,1)]); % pad with zeros for fft
		tmp = fftshift(ifft(atmpi.*conj(tmp)));
		c(:,colaux) = tmp((1:mr)+nfft/2-nr+1);
		ci(col) = i;
		cj(col) = j;
		ci(colaux) = j;
		cj(colaux) = i;
		if ~any(any(imag(a(:,j)))) & real1
			c(:,colaux) = real(c(:,colaux));
		end
		if i~= j
			c(:,col) = conj(c(mr:-1:1,colaux));
		end
	end
end
if nopt == 1	% return normalized by sqrt of each autocorrelation at 0 lag
% do column arithmetic to get correct autocorrelations
	cdiv = ones(mr,1)*sqrt(c(nr,1+(ci-1)*(nc+1)).*c(nr,1+(cj-1)*(nc+1)));
	c = c ./ cdiv;
elseif nopt == 2	% biased result, i.e. divide by nr for each element
	c = c / nr;
elseif nopt == 3	% unbiased result, i.e. divide by nr-abs(lag)
	c = c ./ ([1:nr (nr-1):-1:1]' * ones(1,nsq));
end
if onearray == 1
	c = c(:,1);	% just want the autocorrelation
	[am, an] = size(a);
	if am == 1
		c = c.';
	end
elseif onearray == 2	% produce only cross-correlation
	c = c(:,2);
	if ar == 1
		c = c.';
	end
end
if ~any(any(imag(a)))
	c = real(c);
end

