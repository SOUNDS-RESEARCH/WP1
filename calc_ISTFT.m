function x = calc_ISTFT(X,nfft,noverlap)
%CALC_ISTFT inverse short-time fourier transform using OLA. The ISTFT uses
%a sqrt(hann(nfft)) window.
%
% INPUT:
%   X           : input matrix (mics x bins x frames)
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%
% OUTPUT:
%   x           : output time signal(s)
%
% See also: calc_STFT, ifft.

% Author: Daniel Marquardt & Nico Goessling
% Date: 27.05.2016

if nargin < 3
    noverlap = 2;
end

X = permute(X, [2 3 1]);

% Synthesis window
window  = sqrt(hann(nfft, 'periodic'));

L = size(X,2);
x_tmp = real(ifft([X; conj(X(end-1:-1:2,:,:))], [], 1));

% Apply synthesis window
win_tmp = repmat(window, [1, size(x_tmp, 2), size(x_tmp, 3)]);
x_tmp = x_tmp.*win_tmp;
x_tmp = x_tmp(1:nfft,:,:);

% IFFT per frame
x = zeros((nfft / noverlap)*(L-1) + nfft, size(x_tmp, 3));

% OLA processing
for m = 0:L-1
    x(floor(m * (nfft / noverlap) + 1):floor(m * (nfft / noverlap) + nfft),:) = ...
        squeeze(x_tmp(:,m+1,:)) + ...
        x(floor(m * (nfft / noverlap) + 1):floor(m * (nfft / noverlap) + nfft),:);
end

end