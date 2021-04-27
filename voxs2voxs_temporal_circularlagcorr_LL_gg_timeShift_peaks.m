function voxs2voxs_temporal_circularlagcorr_LL_gg_timeShift_peaks

loc='cluster';
set_parameters;
timeUnit='tr' ;
permN=10000;
lags=-20:20;

for ei=1;%
    
    exp=experiments{ei};
    fs=dir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/isfc_voxind*.mat'  ]);
    fs={fs.name};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/' fs{1}],'r','keptvox','keptT','voxind');
    rz=atanh(r);
    keptvox_target=keptvox;
    [~,tn]=size(rz);
    t_real=(tn-1)/2+1;
    ts_shift=1:tn;
    ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
    tis=  randsample(length(ts_shift),permN,1);
    ts_shift=ts_shift(tis);
    
    peaks=nan(length(fs),length(keptvox_target));
    peakLags=nan(length(fs),length(keptvox_target));
    widths=nan(length(fs),length(keptvox_target));
    peaks_pfwe=nan(length(fs),length(keptvox_target));
    peakLags_pfwe=nan(length(fs),length(keptvox_target));
    widths_pfwe=nan(length(fs),length(keptvox_target));
    peaks_p=nan(length(fs),length(keptvox_target));
    peakLags_p=nan(length(fs),length(keptvox_target));
    widths_p=nan(length(fs),length(keptvox_target));
    rz=[];
    p=[];
    pfwe=[];
    
    for sdi=4475:length(fs);
        try
            load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/perm/' fs{sdi}],'r','keptvox','keptT');
            keptvox_seed=keptvox;
            
            peaks_shift=nan(1,length(keptvox_target),permN);
            for perm=1:permN;
                t_shift=ts_shift(perm);
                [peaks_shift(1,:,perm)]=max(atanh(r(:,t_shift+lags)),[],2);
            end
            
            rz(sdi,:,:)=atanh(r(:,t_real+lags));
            p(sdi,:,:)=mean(permute(repmat(peaks_shift,1,1,1,length(lags)),[1 2 4 3])>rz(sdi,:,:),4);
        catch
            continue
        end
    end
    pfwe=p*length(p(:))/length(lags);
    
    for sdi=1:length(fs);
        for tgi=1:length(keptvox_target);
            temp=squeeze(rz(sdi,tgi,:));
            [pks, locs, ws]=findpeaks(temp,'WidthReference','halfheight');
            locs=locs(pks>0);
            pks=pks(pks>0);
            ws=ws(pks>0);
            if ~isempty(pks);
                [~,loci]=max(pks);
                peakLags(sdi,tgi)=lags(locs(loci));
                peaks(sdi,tgi)=temp(locs(loci));
                widths(sdi,tgi)=ws(loci);
                
                if p(sdi,tgi,locs(loci))<0.05;
                    peakLags_p(sdi,tgi)=lags(locs(loci));
                    peaks_p(sdi,tgi)=temp(locs(loci));
                    widths_p(sdi,tgi)=ws(loci);
                end
                
                if pfwe(sdi,tgi,locs(loci))<0.05;
                    peakLags_pfwe(sdi,tgi)=lags(locs(loci));
                    peaks_pfwe(sdi,tgi)=temp(locs(loci));
                    widths_pfwe(sdi,tgi)=ws(loci);
                end
            end
        end
    end
    
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2voxs//LL_gg/voxs2voxs_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeShift_peaksMax'  ],...
        'rz','lags','keptvox_target','keptvox_seed','keptT','p','pfwe','peakLags','peaks','peaks_pfwe','peakLags_pfwe','widths','widths_pfwe','-7.3');
end



