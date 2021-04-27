%% find the peak nearest to lag 0 instead of the absolute peak
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;

lags_tested={-40:40,-10:10,-15:15, -20:20};
permN=10000;
roi='isc_peak_precuneus';
for ei=1%:11;
    exp=exp_parameters.experiments{ei};
    
    for lagi=4%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/'  roi '.mat'   ],'r','keptvox','keptT');
        rz=atanh(r);
        rzm=nanmean(rz,3);
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/' roi '_timeReversed.mat'   ],'r','keptvox','keptT');
        rz_timeReversed=atanh(r);
        rzm_timeReversed=nanmean(rz_timeReversed,3);
        
        [~,tn,listenerN]=size(rz);
        tmid=(tn-1)/2+1;
        p=[];
        z=[];
        % paired t-test
        
        for vi=1:length(keptvox);
            null=(squeeze(rzm_timeReversed(vi,:)));
            null_m=mean(null);
            null_std=std(null);
            for lagi=1:length(lags);
                r_real=squeeze(rzm(vi,tmid+lags(lagi)));
                [~,p(vi,lagi),~,z(vi,lagi)] = ztest(r_real,null_m, null_std,'tail','right');
            end
        end
        
        [~,~,pfdr]=fdr(p(:),.05);
        pfdr=reshape(pfdr,size(p));
        pfwe=p*length(p(:));
        
        peaks=nan(size(rzm,1),1);
        peakLags=nan(size(rzm,1),1);
        p_peaks=nan(size(rzm,1),1);
        pfdr_peaks=nan(size(rzm,1),1);
        z_peaks=nan(size(rzm,1),1);
        
        for tgi=1:size(rzm,1);
            temp=squeeze(nanmean(rz(tgi,tmid+lags,:),3));
            [pks]=findpeaks(temp);
            pks=pks(pks>0);
            
            if ~isempty(pks);
                pk=max(pks);
                [lagi]=find(temp==pk);
                peakLags(tgi)=lags(lagi);
                peaks(tgi)=pk;
                p_peaks(tgi)=p(tgi,lagi);
                pfdr_peaks(tgi)=pfdr(tgi,lagi);
                z_peaks(tgi)=z(tgi,lagi);
            end
        end
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/'  roi '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks' ],...
            'keptvox','rzm','rz','lags','keptT','p_peaks','p','peakLags','peaks','z_peaks','pfdr','pfdr_peaks');
        
        sig=pfdr_peaks<(.05);
        
        mat=nan(voxn,1);
        mat(keptvox(sig))=peaks(sig);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/' roi '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_fdr_peaks.nii' ]);
        
        mat=nan(voxn,1);
        mat(keptvox(sig))=peakLags(sig);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/' roi '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReverse_fdr_peakLags.nii' ])
        
        sig=p_peaks<(.05/length(p(:)));
        
        mat=nan(voxn,1);
        mat(keptvox(sig))=peaks(sig);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/'  roi '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_fwe_peaks.nii' ]);
        
        mat=nan(voxn,1);
        mat(keptvox(sig))=peakLags(sig);
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/LL_leave1out/'  roi '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_fwe_peakLags.nii' ]);
        
        unique(peakLags(sig))
    end
end


