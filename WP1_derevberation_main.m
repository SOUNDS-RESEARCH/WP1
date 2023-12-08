%WP1 Deliverable
clear all;

%choose this as you desire, currently the folder structure is set such that
%/MATLAB/ is the base directory (delimiter)
str = pwd;
delimiter = 'MATLAB';
parts = strsplit(str, delimiter);
base_path = parts{1};

addpath(genpath(strcat(base_path,delimiter,"\Utilities")));

%This requires to generate the audio data first
load(strcat(base_path,delimiter,"\Audio_Data\audio_data_ICASSP_h_src_1_24.mat"));
load(strcat(base_path,delimiter,"\Audio_Data\audio_data_ICASSP_h_src_25_48.mat"));
load(strcat(base_path,delimiter,"\Audio_Data\audio_data_ICASSP_h_src_dir.mat"));
load(strcat(base_path,delimiter,"\Audio_Data\audio_data_ICASSP_mic_pos.mat"));
load(strcat(base_path,delimiter,"\Audio_Data\audio_data_ICASSP_samp_delay.mat"));
load(strcat(base_path,delimiter,"\Audio_Data\audio_data_ICASSP_src_pos_mat.mat"));
load(strcat(base_path,delimiter,"\Audio_Data\audio_data_ICASSP_target_spk_clean.mat"));

%join together h_src files
h_src = cell(48,3);
h_src(1:24,:) = h_src_1_24;
h_src(25:48,:) = h_src_25_48;


%signal parameters
fs = 16000;

%STFT parameters for WPE
lFrame = 0.064*fs;
L_shift = 0.016*fs;
nb_overlap = lFrame/L_shift;

%compile parameters into STFT_params struct for compact use in functions
STFT_params.lFrame = lFrame;
STFT_params.L_shift = L_shift;
STFT_params.nb_overlap = nb_overlap;
STFT_params.fs = fs;


%STFT parameters for GCC-PHAT
N_gcc = 2048;

%WPE parameters
p = 0.5;
Lg = 8;
tau = 2;
Lg_vec = [8,12,16];
ref_chan = 1;


%Tools for GCC-PHAT estimation
nChan = size(h_src{1},2);
other_mics = 1:nChan;
other_mics(other_mics == ref_chan) = [];

%Parameters for scenario
nSrc_pos = 48;
nRT_60 = 3;

%Initialize cell arrays
perf_reference = cell(nSrc_pos,nRT_60);
WPE_perf_MI = cell(nSrc_pos,nRT_60);
WPE_perf_MD_NINT_ora = cell(nSrc_pos,nRT_60);
WPE_perf_MD_NINT_B2B_ora = cell(nSrc_pos,nRT_60);
WPE_perf_MD_INT_ora = cell(nSrc_pos,nRT_60);
WPE_perf_MD_NINT_estim = cell(nSrc_pos,nRT_60);
WPE_perf_MD_NINT_B2B_estim = cell(nSrc_pos,nRT_60);
WPE_perf_MD_INT_estim = cell(nSrc_pos,nRT_60);



%Metrics booleans
comp_SRMR = false;
comp_CD = true;
comp_LLR = false;
comp_PESQ = true;



%Tool for normalization
length_modified_signal = 176128;




count = 1;
for src_pos_idx = 1:48
    disp(src_pos_idx)
    for rt60_idx = 1:3
        disp(rt60_idx)
        %read audio
        audio_in = fftfilt(h_src{src_pos_idx,rt60_idx},target_spk_clean); %reverberant
        audio_in_dir = fftfilt(h_src_dir{src_pos_idx,rt60_idx},target_spk_clean); %direct component
        %normalize audio (reverberant)
        audio_in_src = audio_in;
        audio_in_src = audio_in_src(1:length_modified_signal,:);
        audio_in_src = audio_in_src./sqrt(var(audio_in_src(:,1)));
        %normalize audio (direct component)
        audio_in_dir_src = audio_in_dir;
        audio_in_dir_src = audio_in_dir_src(1:length_modified_signal,:);
        audio_in_dir_src = audio_in_dir_src./sqrt(var(audio_in_dir_src(:,1)));
        
        %Oracle TDOA calculation
        toa = sqrt(sum((repmat(src_pos_mat(src_pos_idx,:),nChan,1) - mic_pos).^2,2))/340;
        tdoa = toa - toa(ref_chan);
        TDOA = -1*tdoa*fs;

        %Estimated TDOA calculation
        TDOA_temp = zeros(nChan,343);
        for i = 1:length(other_mics)
            other_mic = other_mics(i);
            TDOA_temp(other_mic,:) = (gcc_batch(audio_in_src(:,[ref_chan,other_mic]),N_gcc,1,1,N_gcc/2,N_gcc*2,10)).';
        end
        TDOA_estim = mode(TDOA_temp,2);
        
        %Perform STFT
        [X,~] = calc_STFT(audio_in_src,fs,lFrame,nb_overlap);
        X = permute(X,[1 3 2]);
        
        [~,nFrames,nFreq] = size(X); %Matrix X has dimension nChan x nFrames x nFreq
        
        %Compute reference performance
        perf_reference{src_pos_idx,rt60_idx} = comp_perf(audio_in_src(:,ref_chan),audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);

        %% baseline fixed (microphone-independent) prediction delay

        d_WPE_MI = WPE_MDPD(X,p,tau,Lg_vec(rt60_idx),ref_chan,'MI',[],STFT_params);
    
        audio_out_MI = calc_ISTFT(permute(d_WPE_MI,[3 2 1]),lFrame,nb_overlap);

        WPE_perf_MI{src_pos_idx,rt60_idx} = comp_perf(audio_out_MI,audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);


        %% microphone-dependent prediction delay (oracle TDOA)
        
        d_WPE_MD_NINT_ora = WPE_MDPD(X,p,tau,Lg_vec(rt60_idx),ref_chan,'MD-NINT',TDOA,STFT_params);

        audio_out_MD_NINT_ora = calc_ISTFT(permute(d_WPE_MD_NINT_ora,[3 2 1]),lFrame,nb_overlap);

        WPE_perf_MD_NINT_ora{src_pos_idx,rt60_idx} = comp_perf(audio_out_MD_NINT_ora,audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);


        
        
        d_WPE_MD_NINT_B2B_ora = WPE_MDPD(X,p,tau,Lg_vec(rt60_idx),ref_chan,'MD-NINT-B2B',TDOA,STFT_params);
        
        audio_out_MD_NINT_B2B_ora = calc_ISTFT(permute(d_WPE_MD_NINT_B2B_ora,[3 2 1]),lFrame,nb_overlap);

        WPE_perf_MD_NINT_B2B_ora{src_pos_idx,rt60_idx} = comp_perf(audio_out_MD_NINT_B2B_ora,audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);
        
        
        
        
        
        d_WPE_MD_INT_ora = WPE_MDPD(X,p,tau,Lg_vec(rt60_idx),ref_chan,'MD-INT',TDOA,STFT_params);
        
        audio_out_MD_INT_ora = calc_ISTFT(permute(d_WPE_MD_INT_ora,[3 2 1]),lFrame,nb_overlap);

        WPE_perf_MD_INT_ora{src_pos_idx,rt60_idx} = comp_perf(audio_out_MD_INT_ora,audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);

    
    
        %% microphone-dependent prediction delay (estimated TDOA)
   
        
        d_WPE_MD_NINT_estim = WPE_MDPD(X,p,tau,Lg_vec(rt60_idx),ref_chan,'MD-NINT',TDOA_estim,STFT_params);

        audio_out_MD_NINT_estim = calc_ISTFT(permute(d_WPE_MD_NINT_estim,[3 2 1]),lFrame,nb_overlap);

        WPE_perf_MD_NINT_estim{src_pos_idx,rt60_idx} = comp_perf(audio_out_MD_NINT_estim,audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);


        
        
        d_WPE_MD_NINT_B2B_estim = WPE_MDPD(X,p,tau,Lg_vec(rt60_idx),ref_chan,'MD-NINT-B2B',TDOA_estim,STFT_params);
        
        audio_out_MD_NINT_B2B_estim = calc_ISTFT(permute(d_WPE_MD_NINT_B2B_estim,[3 2 1]),lFrame,nb_overlap);

        WPE_perf_MD_NINT_B2B_estim{src_pos_idx,rt60_idx} = comp_perf(audio_out_MD_NINT_B2B_estim,audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);
        
        
        
        
        
        d_WPE_MD_INT_estim = WPE_MDPD(X,p,tau,Lg_vec(rt60_idx),ref_chan,'MD-INT',TDOA_estim,STFT_params);
        
        audio_out_MD_INT_estim = calc_ISTFT(permute(d_WPE_MD_INT_estim,[3 2 1]),lFrame,nb_overlap);

        WPE_perf_MD_INT_estim{src_pos_idx,rt60_idx} = comp_perf(audio_out_MD_INT_estim,audio_in_dir_src(:,ref_chan),fs,comp_SRMR,comp_PESQ,comp_LLR,comp_CD);
    
    end
end

