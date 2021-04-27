function voxs2voxs_temporal_circularlagcorr_LL_gg_timeShift_peaksNearest

loc='cluster';
set_parameters;
timeUnit='tr' ;
permN=10000;
lags=-15:15;

for ei=1;%
    
    exp=experiments{ei};
    fs=dir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/isfc_voxind*.mat'  ]);
    fs={fs.name};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/' fs{1}],'r','keptvox','keptT','voxind');
    rz=atanh(r); clear r
    
    [~,tn]=size(rz);
    t_real=(tn-1)/2+1;
    ts_shift=1:tn;
    ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
    tis=  randsample(length(ts_shift),permN,1);
    ts_shift=ts_shift(tis);
    
    peaks=nan(length(fs),length(keptvox));
    peakLags=nan(length(fs),length(keptvox));
    p=nan(length(fs),length(keptvox),length(lags));
    pfwe=nan(length(fs),length(keptvox),length(lags));
    peaks_pfwe=nan(length(keptvox),length(keptvox));
    peakLags_pfwe=nan(length(keptvox),length(keptvox));
    peaks_p=nan(length(keptvox),length(keptvox));
    peakLags_p=nan(length(keptvox),length(keptvox));
    rz_real=nan(length(fs),length(keptvox),length(lags));
    for sdi=1:length(fs);
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/' fs{sdi}],'r','keptvox','keptT');
        rz=atanh(r); clear r
        keptvox=keptvox;
        
        peaks_shift=[];
        for perm=1:permN;
            t_shift=ts_shift(perm);
            [peaks_shift(:,perm)]=max(rz(:,t_shift+lags),[],2);
        end
        
        rz=rz(:,t_real+lags);
        rz_real(sdi,:,:)=rz;
        p(sdi,:,:)=mean(permute(repmat(peaks_shift,1,1,length(lags)),[1 3 2])>rz,3);
        pfwe(sdi,:,:)=p(sdi,:,:)*length(p(:))/length(lags);
        
        for ni=1:length(keptvox);
            temp=squeeze(rz(ni,:));
            [pks, locs]=findpeaks(temp);
            locs=locs(pks>0);
            pks=pks(pks>0);
            if ~isempty(pks);
                [~,loci]=(min(abs(locs-find(lags==0))));
                peakLags(sdi,ni)=lags(locs(loci));
                peaks(sdi,ni)=temp(locs(loci));
                
                if p(sdi,ni,locs(loci))<0.05;
                    peakLags_p(sdi,ni)=lags(locs(loci));
                    peaks_p(sdi,ni)=temp(locs(loci));
                end
                
                if pfwe(sdi,ni,locs(loci))<0.05;
                    peakLags_pfwe(sdi,ni)=lags(locs(loci));
                    peaks_pfwe(sdi,ni)=temp(locs(loci));
                end
            end
        end
    end
    
    rz=rz_real;
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/voxs2voxs_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksNearest'  ],...
       'lags','keptvox','keptvox','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe','-v7.3')
end


