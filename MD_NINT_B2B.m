function[X_postproc] = MD_NINT_B2B(TDOA,X,tau,STFT_params)


[nChan,nFrames,~] = size(X); %dimensions of X should be nChan x nFrames x nFreq


%The way this script is structured is:
%   1. Compute non-integer frame delay time-domain TDOA compensation filter 
%   2. Translate time-domain TDOA compensation filter into crossband STFT
%   filter u_m
%   3. Filter input signal with crossband filter u_m and integer delay
%   tau_int


%taking values from STFT params struct
L_shift = STFT_params.L_shift;
lFrame = STFT_params.lFrame;
nb_overlap = STFT_params.nb_overlap;
fs = STFT_params.fs;


%some parameter definitions
filter_length = 2*L_shift + 1;
nFreq_filter = lFrame/2 + 1;
lagr_poly_deg = 8;


    
%TDOA components
delta_frame = round(TDOA/L_shift);
delta_samp = floor(TDOA - (delta_frame)*L_shift);
delta_frac = TDOA - delta_frame*L_shift - delta_samp;



%1. Compute non-integer frame delay time-domain TDOA compensation filter

%1.1 Setting up sinc function (fractional compensation filter)
FIR_length = lagr_poly_deg+1;  %Sinc FIR Length
idxWindow = (-floor((FIR_length-1)/2):floor(FIR_length/2))'; %Corresponding indexes (given acausal sinc)
midpoint = ceil(FIR_length/2); %FIR filter midpoint (origin)

% Define and solve Lagrange interpolation equations
V = idxWindow.^(0:lagr_poly_deg); % Vandermonde structure

C = delta_frac.^(0:lagr_poly_deg);

frac_delay_filt = C/V;  % Solve for the coefficients

%1.2 Setting up integer delay filter

%causal shifts
causal_shift_integer = L_shift + 1;
causal_shift_frac_filter = midpoint;


samp_delay_filt = zeros(filter_length,nChan);

for i = 1:nChan
    samp_delay_filt(delta_samp(i) + causal_shift_integer,i) = 1; %make it such that positive is delay
end
%1.3 Combining integer delay filter and fractional delay filter to generate
%remaining time-domain TDOA compensation filter
nonintegerframe_TDOA_comp_filt = zeros(filter_length + FIR_length - 1,nChan);
for i = 1:nChan
    nonintegerframe_TDOA_comp_filt(:,i) = conv(samp_delay_filt(:,i),frac_delay_filt(i,:).');
end


%window function used (modify accordingly)
window  = sqrt(hann(lFrame, 'periodic'));

%2.1 compute phi
phi_term = comp_phi(lFrame,nFreq_filter,window);

causal_shift_window_sum_m_term = lFrame;

causal_shift_total = causal_shift_integer + causal_shift_window_sum_m_term + causal_shift_frac_filter;

%2.2 compute u_m
u_m = comp_u_m(nonintegerframe_TDOA_comp_filt,phi_term,causal_shift_total,L_shift,nChan,nFreq_filter);

%3. filtering step
length_u_m = size(u_m,3);
filter_offset = (length_u_m-1)/2;

tau_int = tau + delta_frame;


X_postproc = zeros(nChan,nFrames,nFreq_filter);

for i = 1:nChan
    X_temp = zeros(nFrames,nFreq_filter);  
    for k = 1:nFreq_filter
        u_m_temp = squeeze(u_m(k,:,i));
        X_temp_conv = squeeze(X(i,:,k));
        X_temp_temp = conv(u_m_temp,X_temp_conv).';
        X_temp(:,k) = X_temp_temp(1+filter_offset:filter_offset + nFrames,:);
    end
    X_postproc_temp = X_temp;
    if tau_int(i) >= 0
        X_postproc(i,:,:) = [zeros(tau_int(i),nFreq_filter);X_postproc_temp(1:end-tau_int(i),:)];
    elseif tau_int(i) < 0
        X_postproc(i,:,:) = [X_postproc_temp(1+tau_int(i):end,:);zeros(tau_int(i),nFreq_filter)];
    end
end
%optional normalization (used to match normalization done on speech signals
%earlier)
X_postproc = normalize(X_postproc,lFrame,nb_overlap,fs);

end

function phi_term = comp_phi(lFrame,nFreq_filter,window)

    e_term = exp(1i*2*pi*((0:(lFrame/2)).'/lFrame)*(-(lFrame-1):(lFrame-1)));
    
    
    
    window_sum_m_term = zeros(nFreq_filter,2*lFrame-1);

    k_dash_offset_range = 0;

        
    for k = 1:nFreq_filter
        k_dash_offset_range_temp = k_dash_offset_range;
    
        k_dash_offset_range_temp(k+k_dash_offset_range_temp < 1 | k+k_dash_offset_range_temp > nFreq_filter) = [];
        k_dash_indexs = k + k_dash_offset_range_temp;
    
        window_sum_m_term(k,:) = conv2(window.*exp(-1*1i*((2*pi)/lFrame)*(0:lFrame-1).'*k_dash_offset_range_temp),window).';
    end
    
    phi_term = zeros(size(window_sum_m_term));
    
    for k = 1:nFreq_filter
        phi_term(k,:) = window_sum_m_term(k,:).*e_term(k,:);
    end
end

function u_m = comp_u_m(td_filter,phi_term,causal_shift_total,L_shift,nChan,nFreq_filter)

    p_0th_index = causal_shift_total - 1 - 1;
    
    
    n_taps = round(2*(p_0th_index/L_shift - 1) + 1);
    tap_index_offset = (-(n_taps - 1)/2:(n_taps - 1)/2)*L_shift;
    tap_indexs = p_0th_index + tap_index_offset;

    u_m = zeros(nFreq_filter,n_taps,nChan);
    for i = 1:nChan
        td_filter_reshape = reshape(td_filter(:,i),[1,size(td_filter,1)]);
        h_term = conv2(phi_term,td_filter_reshape);   
    
        u_m_temp = h_term(:,tap_indexs);

        u_m(:,:,i) = u_m_temp;
    end
end


function X_normalized = normalize(X,lFrame,nb_overlap,fs)
    X = permute(X,[1 3 2]);

    audio_out_X = calc_ISTFT(X,lFrame,nb_overlap);
    
    audio_out_X = audio_out_X./sqrt(var(audio_out_X(:,1)));
    [X_normalized,~] = calc_STFT(audio_out_X,fs,lFrame,nb_overlap);
    X_normalized = permute(X_normalized,[1 3 2]);
end