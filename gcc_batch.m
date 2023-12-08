function [estdel,ccw] = gcc_batch(y,N,weighting,wintype,overlap,L,intfactor)
%
% Compute delay between 2 signals using Generalised Crosscorrelation
% Using zero-padding in frames (circular convolution effects)
%
% Reference: C.H. Knapp and G.C. Carter, "The Generalized Correlation Method for Estimation of Time Delay," 
%            IEEE Trans. Acoust., Speech and Signal Proc., vol. 24, no. 4, Aug. 1976, pp. 320-327.
%
% OUTPUT    estdel         Estimated delay (number of samples)
%           ccw            Weighted time-domain cross-correlation  
%
% INPUTS    y              signal
%           N              framesize
%           lambda         averaging factor for spectra (lambda = 0, no averaging)
%           weighting      Weighting function (optional)
%                            0: No weighting
%                            1: PHAT, Phase Transform (default)
%                            2: ML, Maximum Likelihood
%                            3: SCOT, Smoothed Coherence Transform
%                            4: Roth
%                            5: Eckart
%           wintype     Windowing function for computing psd/coherence (optional)
%                            0: Rectangular window (default)
%                            1: Hanning window
%                            2: Hamming window
%           overlap        number of overlap samples (optional, default N/2)
%           L              size of FFT (optional, default 2^(nextpow2(N)+1))
%           intfactor      Interpolation factor for GCC function (optional, default 1)
%                            
  


totallength = size(y,1);  
  
if nargin < 7,
  intfactor = 1;
  if nargin < 6,
    L = 2^(nextpow2(N)+1);
    if nargin < 5,
      overlap = N/2;
      if nargin < 4,
        wintype = 0;
        if nargin < 3,
          weighting = 1;
        end
      end
    end
  end
end

% Windowing
if wintype == 0,
  win = ones(N,1);
elseif wintype == 1,
  win = hanning(N);
elseif wintype == 2,
  win = hamming(N);
end

NrFrames = ceil((totallength-N)/overlap);
framestart = 1;

C1 = zeros(L,1);
C2 = zeros(L,1);
CC = zeros(L,1);
coh = zeros(L,1);
coh2 = zeros(L,1);
ccw = zeros(L*intfactor,1);
CCw = zeros(L,1);
%estdel = zeros(1);

% Initialisation using first frame
frame1 = y(framestart:framestart+N-1,1).*win;
frame2 = y(framestart:framestart+N-1,2).*win;
F1(:,1) = fft(frame1,L);
F2(:,1) = fft(frame2,L);
C1(:,1) = abs(F1.*conj(F1)); % PSD1 (real)
C2(:,1) = abs(F2.*conj(F2)); % PSD2 (real)
CC(:,1) = F1.*conj(F2); % Cross-correlation
framestart = framestart+overlap;

for i=2:NrFrames,
  
  frame1 = y(framestart:framestart+N-1,1).*win;
  frame2 = y(framestart:framestart+N-1,2).*win;

  % Frequency-domain correlation
  F1(:,i) = fft(frame1,L);
  F2(:,i) = fft(frame2,L);
    
  framestart = framestart+overlap;

end

C1 = (1/NrFrames)*sum(F1.*conj(F1),2);
C2 = (1/NrFrames)*sum(F2.*conj(F2),2);
CC = (1/NrFrames)*sum(F1.*conj(F2),2);
coh = CC./sqrt(C1.*C2);



  % Weighting
  if weighting == 0,
    weight = ones(L,1);
  elseif weighting == 1,
    weight = 1./abs(CC);
  elseif weighting == 2,
    coh2 = abs(coh).^2; % Magnitude-squared coherence
    weight = coh2./(abs(CC).*(1-coh2));
  elseif weighting == 3,
    weight = 1./sqrt(C1.*C2);
  elseif weighting == 4,
    weight = 1./C1;
  elseif weighting == 5,
    weight = abs(CC)./((C1-abs(CC)).*(C2-abs(CC)));
  end

  CCw = CC.*weight;

  % Estimation of delay

  tmp = fftshift(real(ifft(CCw)));

  if intfactor == 1;
    ccw = tmp;
    [m,index] = max(ccw);
    estdel = index-L/2-1;
  else
    ccw = interp(tmp,intfactor);
    [m,index] = max(ccw);
    estdel = (index-1)/intfactor-L/2;
  end



end
