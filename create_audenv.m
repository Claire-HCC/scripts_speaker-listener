%% stim power vs BOLD
clear
set_parameters;

resamp_factor = 1;
cropped='';%'_uncropped',"";
role='listener';

for ei = 6;
    exp=exp_parameters.experiments{ei};
     stimname =ls([expdir exp '/sound/*wav']);
    stimname=[exp '_' role cropped];
    % stimname=ls([expdir exp '/sound/*wav']);
    % stimname=strrep( stimname,'.wav','');
    
    [w, audfs] = audioread([expdir exp '/sound/' stimname '.wav']);
    
    w0=w;
    %audiowrite([expdir experiments{ei} '/sound/' experiments{ei} '.wav'],w0,audfs)
    
    w=w0(:,1);
    w = resample(w, 1, resamp_factor);  %downsample for easier manipulation
    
    audfs = audfs/resamp_factor;
    audtime = 0 + ( 0 : (length(w)-1) ) / audfs;
    
    q = w;
    stim_mod = abs(hilbert(q));
    stim_time = audtime;
    
    boldcrop = stim_time > 0; % wav_crop_start(ei);
    audpow = stim_mod(boldcrop);
    audpow_slow = resample(audpow,1,(audfs/resamp_factor)*exp_parameters.tr(ei));
    
    tshift = 0; %
    audpow_slow = round(1000+100*zscore(audpow_slow));
    audpow_slow_shifted_100 = round(1000+100*zscore([audpow_slow(1:tshift); audpow_slow(1:end-tshift)]));
    audenv = audpow_slow_shifted_100;
    aud=audenv;
    
    save( [expdir exp '/sound/'  stimname '_audenv'],'aud');
    %   save([expdir experiments{ei} '/sound/'  experiments{ei} '_' role '_'  cropped '_audenv'],'aud');
    
    xBF.dt=exp_parameters.tr(ei); % TR
    xBF.name='hrf'% (with time and dispersion derivatives)';
    bf = spm_get_bf(xBF)
    
    audenv=audenv-mean(audenv);
    aud=conv(audenv,bf.bf);
    save( [expdir exp '/sound/'  stimname '_audhrf'],'aud');
    %  save( [expdir experiments{ei} '/sound/'  experiments{ei} '_' role '_' cropped '_audhrf'],'aud');
    
end



