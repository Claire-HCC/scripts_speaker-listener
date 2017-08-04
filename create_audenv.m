%% stim power vs BOLD
clear
set_parameters;

resamp_factor = 1;

for ei = 1:2;
  %  stimname =ls([expdir experiments{ei} '/sound/*wav']);
   stimname=experiments{ei};
    stimname=strrep([expdir experiments{ei} '/sound/' stimname],'.wav','');
  
    [w, audfs] = audioread([stimname '.wav']);
    
    w0=w(wav_crop_start*audfs:end);
    audiowrite([expdir experiments{ei} '/sound/' experiments{ei} '.wav'],w0,audfs)
    
    w=w0(:,1);
    w = resample(w, 1, resamp_factor);  %downsample for easier manipulation
    
    audfs = audfs/resamp_factor;
    audtime = 0 + ( 0 : (length(w)-1) ) / audfs;
    
    q = w;
    stim_mod = abs(hilbert(q));
    stim_time = audtime;
    
    boldcrop = stim_time > 0; % wav_crop_start(ei);
    audpow = stim_mod(boldcrop);
    audpow_slow = resample(audpow,1,(audfs/resamp_factor)*tr(ei));
    
    tshift = 0; %
    audpow_slow = round(1000+100*zscore(audpow_slow));
    audpow_slow_shifted_100 = round(1000+100*zscore([audpow_slow(1:tshift); audpow_slow(1:end-tshift)]));
    audenv = audpow_slow_shifted_100;
    aud=audenv;
    
    save([ experiments{ei}  '_audenv'],'aud');
    
    xBF.dt=tr(ei); % TR
    xBF.name='hrf'% (with time and dispersion derivatives)';
    bf = spm_get_bf(xBF)
    
    audenv=audenv-mean(audenv);
    aud=conv(audenv,bf.bf);
    save([ experiments{ei}  '_audhrf'],'aud');
    
end



