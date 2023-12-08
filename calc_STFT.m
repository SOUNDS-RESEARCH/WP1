function [X,f] = calc_STFT(x,fs,nfft,noverlap)
%CALC_STFT short-time fourier transform using OLA. The STFT uses a
%sqrt(hann(nfft)) window.
%
% INPUT:
%   x           : input time signal(s) (samples x channels)
%   fs          : sampling rate
%   nfft        : FFT size
%   noverlap    : frame overlap; default: 2 (50%)
%
% OUTPUT:
%   X           : STFT matrix (channels x bins x frames)
%   f           : frequency vector for bins
%
%   See also: calc_ISTFT, fft.

% Author: Daniel Marquardt & Nico Goessling
% Date: 27.05.2016

if nargin < 4
    noverlap = 2;
end

% synthesis window
window  = sqrt(hann(nfft, 'periodic'));

% use only half FFT spectrum
N_half = nfft / 2 + 1;

% get frequency vector
f = 0:(fs / 2) / (N_half - 1):fs / 2;

% init
L = floor((length(x) - nfft + (nfft / noverlap)) / (nfft / noverlap));
X = zeros(N_half, L, size(x,2));

% OLA processing
for l = 0:L-1 % Frame index
    x_frame = x(floor(l*(nfft / noverlap) + 1):floor(l*(nfft / noverlap) + nfft),:);
    x_windowed = x_frame.*repmat(window, 1, size(x_frame,2));
    X_frame =  fft(x_windowed,[],1);
    X(:,l+1,:) = X_frame(1 : floor(size(X_frame,1) / 2) + 1, :);
end

X = permute(X, [3 1 2]);

end