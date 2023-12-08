function[X_postproc] = MD_INT(TDOA,X,tau,STFT_params)

L_shift = STFT_params.L_shift;

delta_frame = round(TDOA/L_shift);

tau_int = tau + delta_frame;

[nChan,nFrames,nFreq] = size(X);
X_postproc = zeros(nChan,nFrames,nFreq);

for i = 1:nChan
    if tau_int(i) >= 0
        X_postproc(i,:,:) = [zeros(tau_int(i),nFreq);squeeze(X(i,1:end-tau_int(i),:))];
    elseif tau_int(i) < 0
        X_postproc(i,:,:) = [squeeze(X(i,1+tau_int(i):end-tau_int(i),:));zeros(tau_int(i),nFreq)];
    end
end





