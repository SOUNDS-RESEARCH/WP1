function[d_out] = WPE_MDPD(X,p,tau,Lg,ref_chan,mode,TDOA,STFT_params)

%dimensions of X should be nChan x nFrames x nFreq
if (size(X,1) > size(X,2)) || (size(X,1) > size(X,3))
    warning("Number of channels is greater than number of STFT frames or number of Frequency bins, check dimensions")
end

%Additional WPE parameters
maxIt = 10;
reg_factor = 1e-8; 
nChan = size(X,1);
nFrames = size(X,2);
nFreq = size(X,3);

%compute prediction signal matrix Y depending on choice of mode
if strcmp(mode,'MD-NINT')
    Y = MD_NINT(TDOA,X,tau,STFT_params);
elseif strcmp(mode,'MD-NINT-B2B')
    Y = MD_NINT_B2B(TDOA,X,tau,STFT_params);
elseif strcmp(mode,'MD-INT')
    Y = MD_INT(TDOA,X,tau,STFT_params);
elseif strcmp(mode,'MI')
    if tau >= 0
        Y = cat(2,zeros(nChan,tau,nFreq),X(:,1:end-tau,:));
    elseif tau < 0
        Y = cat(2,X(:,1+tau,:),zeros(nChan,tau,nFreq));
    end
else
    error('unrecognized MDPD mode selected')
end


%Due to page-base processing, the first (DC) and last frequency bin are
%removed for processing and added back afterwards
X = X(:,:,2:end-1);
Y = Y(:,:,2:end-1);

%update nFreq
nFreq = size(X,3);

%set reference signal
x_ref = X(ref_chan,:,:);
x_ref = permute(x_ref,[2 1 3]);


%% Build convolutional matrix
R_convmtx = zeros(nChan,Lg,nFrames);
X_tau = zeros(nChan*(Lg),nFrames,nFreq);
for k = 1:size(X,3)
    %build convolutional STFT signal
    for m = 1:nChan
        convmtx_temp = convmtx(Y(m,:,k),Lg);
        R_convmtx(m,:,:) = convmtx_temp(:,1:nFrames);
    end
    X_tau(:,:,k) = reshape(R_convmtx(:),nChan*size(R_convmtx,2),[]);
end

X_tau = permute(X_tau,[2 1 3]);

%% Loop over iterations optimizing the variance coefficients
for iter = 1:maxIt
    if iter ~= 1
        d_temp = d;
        d_temp = permute(d_temp, [2,4,1,3]);
    else
        d_temp = x_ref;
        d_temp = permute(d_temp, [2,4,1,3]);
    end
    
    %compute WPE weights
    WPE_weights = permute((abs(pagemtimes(d_temp, 'ctranspose', d_temp, 'none')) + reg_factor).^(p/2-1), [3,1,4,2]);
    %compute weighted convolution matrix
    weighted_X_tau = WPE_weights.* X_tau;
    %compute covariance and cross covariance matrices
    WPE_CovMat = pagemtimes(weighted_X_tau, 'ctranspose', X_tau, 'none');
    WPE_CrossCovMat = pagemtimes(weighted_X_tau, 'ctranspose', x_ref,'none');
    %perform regularization based on trace
    reg_trace_param = permute(arrayfun(@(subband) trace(WPE_CovMat(:,:,subband)), 1:nFreq), [3,1,2]);
    reg_trace_param(reg_trace_param == 0) = 1;
    WPE_CovMat = WPE_CovMat + reg_factor .* reg_trace_param .* repmat(eye(size(WPE_CovMat, [1,2])), [1,1,nFreq]);

    %compute reverb filter
    reverb_filter = pagemldivide(WPE_CovMat, WPE_CrossCovMat); %requires MATLAB R2022b
    %compute reverb signal
    reverb_stft_signal = pagemtimes(X_tau, 'none', reverb_filter, 'none');
    %compute dereverberated output
    d = x_ref - reverb_stft_signal;

end
%re-insert DC and last subband
d_out = cat(2,squeeze(X(1,:,1)).',squeeze(d),squeeze(X(1,:,end)).');

end
