function voxs2voxs_temporal_circularlagcorr_LL_leave1out_timeShift_peaks

% loc='cluster';
set_parameters;
timeUnit='tr' ;
permN=10000;
lags=-15:15;
for ei=12;%
    
    exp=experiments{ei};
    fs=dir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_leave1out/perm/isfc_voxind*.mat'  ]);
    fs={fs.name};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_leave1out/perm/' fs{1}],'rzm','keptvox','keptT','voxind');
    keptvox_target=keptvox;
    [~,tn]=size(rzm);
    t_real=(tn-1)/2+1;
    ts_shift=1:tn;
    ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
    tis=  randsample(length(ts_shift),permN,1);
    ts_shift=ts_shift(tis);
    
    peaks=nan(length(fs),length(keptvox_target));
    peakLags=nan(length(fs),length(keptvox_target));
    peaks_shift=nan(length(fs),length(keptvox_target),permN);
    
    rzm_real=nan(length(fs),length(keptvox_target),length(lags));
    
    for sdi=1:length(fs);
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_leave1out/perm/' fs{sdi}],'rzm','keptvox','keptT');
        keptvox_seed(sdi,1)=voxind;
        [peaks(sdi,:),lagi]=max(rzm(:,t_real+lags),[],2);
        peakLags(sdi,:)=lags(lagi);
        
        for perm=1:permN;
            t_shift=ts_shift(perm);
            [peaks_shift(sdi,:,perm)]=max(rzm(:,t_shift+lags),[],2);
        end
        
        p=mean(peaks_shift>peaks,3);
        pfwe=p*(size(peaks,1)*size(peaks,2));
        [~,~,pfdr]=fdr(p(:));
        pfdr=reshape(pfdr,size(p));
        
        peakLags_pfwe=peakLags;
        peakLags_pfwe(pfwe>0.05)=NaN;
        peaks_pfwe=peaks;
        peaks_pfwe(pfwe>0.05)=NaN;
        
        peakLags_pfdr=peakLags;
        peakLags_pfdr(pfdr>0.05)=NaN;
        peaks_pfdr=peaks;
        peaks_pfdr(pfdr>0.05)=NaN;
             rzm_real(sdi,:,:)=rzm(:,t_real+lags);
    end
rzm=rzm_real;
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_leave1out/voxs2voxs_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaks'  ],...
            'rzm','lags','keptvox_target','keptvox_seed','keptT','p','pfdr','pfwe','peakLags','peaks','peaks_pfdr','peakLags_pfdr','peaks_pfwe','peakLags_pfwe')
end


