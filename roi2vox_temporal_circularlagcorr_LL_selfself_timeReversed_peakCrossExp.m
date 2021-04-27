% find the peak nearest to lag 0 instead of the absolute peak
clear all

loc='mypc';
set_parameters;

timeUnit='tr' ;
froidir='isc_peak';
lags=-15:15;

rzms=[];
rzms_timeReversed=[];
rois={'Angular_L','precuneus','HG_L'};
perc=0.30;

eis=[1:8];


for ri=1:length(rois);
    roi=rois{ri};
    rzms=nan(voxn,201,length(eis));
    rzms_timeReversed=nan(voxn,201,length(eis));
    
    for ei=1:length(eis);
        exp=exp_parameters.experiments{ei};
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/'  roi '.mat'   ],'r','keptvox','keptT');
        rzm=nanmean(atanh(r),3);
        
        load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/'  froidir '/LL_selfself/' roi '_timeReversed.mat'   ],'r','keptT');
        rzm_timeReversed=nanmean(atanh(r),3);
        
        [~,tn]=size(rzm);
        tmid=(tn-1)/2+1;
        
        rzms(keptvox,:,ei)=rzm(:,tmid+[-100:100]);
        rzms_timeReversed(keptvox,:,ei)=rzm_timeReversed(:,tmid+[-100:100]);
    end
    keptvox=find(sum(isnan(rzms(:,101,:)),3)==0);
    rzms=rzms(keptvox,:,:);
    rzms_timeReversed=rzms_timeReversed(keptvox,:,:);
    rzmsm=nanmean(rzms,3);
    rzmsm_timeReversed=nanmean(rzms_timeReversed,3);
    
    % paired t-test
    tmid=(size(rzms,2)-1)/2+1;
    p_time=[];
    z_time=[];
    % paired t-test
    for vi=1:length(keptvox);
        
        null=(squeeze(rzmsm_timeReversed(vi,:)));
        null_m=mean(null);
        null_std=std(null);
        for lagi=1:length(lags);
            r_real=squeeze(rzmsm(vi,tmid+lags(lagi)));
            [~,p_time(vi,lagi),~,z_time(vi,lagi)] = ztest(r_real,null_m, null_std,'tail','right');
        end
        
    end
    [~,~,pfdr_time]=fdr(p_time(:),.05);
    pfdr_time=reshape(pfdr_time,size(p_time));
    
    p_exps=[];
    t_exps=[];
    peakLags=[];
    for vi=1:length(keptvox);
        for lagi=1:length(lags);
            r_real=squeeze(rzms(vi,tmid+lags(lagi),:))';
            [~,p_exps(vi,lagi),~,stats] = ttest(r_real,0,'tail','right');
            t_exps(vi,lagi)=stats.tstat;
        end
    end
    
    [~,~,pfdr_exps]=fdr(p_exps(:),.05);
    pfdr_exps=reshape(pfdr_exps,size(p_exps));
    
    peaks=nan(length(keptvox),1);
    peakLags=nan(length(keptvox),1);
    p_time_peaks=nan(length(keptvox),1);
    pfdr_time_peaks=nan(length(keptvox),1);
    p_exps_peaks=nan(length(keptvox),1);
    pfdr_exps_peaks=nan(length(keptvox),1);
    t_peaks=nan(length(keptvox),1);
    
    npeaks=nan(length(keptvox),1);
    npeakLags=nan(length(keptvox),1);
    % p_time_npeaks=nan(length(keptvox),1);
    % pfdr_time_npeaks=nan(length(keptvox),1);
    % p_exps_npeaks=nan(length(keptvox),1);
    % pfdr_exps_npeaks=nan(length(keptvox),1);
    % t_npeaks=nan(length(keptvox),1);
    
    for vi=1:length(keptvox);
        
        temp=squeeze(rzmsm(vi,tmid+lags));
        [pks]=findpeaks(temp);
        pks=pks(pks>0);
        
        if ~isempty(pks);
            pk=max(pks);
            [lagi]=find(temp==pk);
            peakLags(vi,1)=lags(lagi);
            peaks(vi,1)=pk;
            p_time_peaks(vi,1)=p_time(vi,lagi);
            pfdr_time_peaks(vi,1)=pfdr_time(vi,lagi);
            p_exps_peaks(vi,1)=p_exps(vi,lagi);
            pfdr_exps_peaks(vi,1)=pfdr_exps(vi,lagi);
            %    t_peaks(vi,1)=t(vi,lagi);
        end
        
        temp=-squeeze(rzmsm(vi,tmid+lags));
        [pks]=findpeaks(temp);
        pks=pks(pks>0);
        if ~isempty(pks);
            pk=-max(pks);
            [lagi]=find(temp==max(pks));
            npeakLags(vi,1)=lags(lagi);
            npeaks(vi,1)=pk;
            %         p_time_npeaks(vi,1)=p_time(vi,lagi);
            %         pfdr_time_npeaks(vi,1)=pfdr_time(vi,lagi);
            %         p_exps_npeaks(vi,1)=p_exps(vi,lagi);
            %         pfdr_exps_npeaks(vi,1)=pfdr_exps(vi,lagi);
            %      t_npeaks(vi,1)=t(vi,lagi);
        end
    end
    
    
    exp='crossExps';
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir  '/LL_selfself/' ]);
    rzm=rzmsm;
    pfdr_peaks=pfdr_exps_peaks;
    pfdr_peaks(:,2)=pfdr_time_peaks;
    pfdr_peaks=max(pfdr_peaks,[],2);
    % pfdr_npeaks=pfdr_exps_npeaks;
    % pfdr_npeaks(:,2)=pfdr_time_npeaks;
    % pfdr_npeaks=max(pfdr_npeaks,[],2);
    
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir  '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReverse_fdr_peakLags.mat' ],...
        '','rzms','rzms_timeReversed','rzmsm','rzmsm_timeReversed','lags','peakLags','peaks','rzm','pfdr_peaks','pfdr_exps_peaks','pfdr_time_peaks','tmid','npeaks');
end

exp='crossExps';
for ri=1:length(rois);
    roi=rois{ri};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir  '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReverse_fdr_peakLags.mat' ],...
        '','rzms','rzms_timeReversed','rzmsm','rzmsm_timeReversed','lags','peakLags','peaks','rzm','pfdr_peaks','pfdr_exps_peaks','pfdr_time_peaks','tmid','npeaks');
    
    sig=(pfdr_peaks<.05 & (abs(npeaks)<peaks | isnan(npeaks) ));
    unique(peakLags(sig))
    
    % mat=zeros(voxn,1);
    % mat(keptvox)=sig;
    % nii=mat2nii(mat);
    % save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_pfdr_mask.nii' ]);
    %
    %
    % mat=nan(voxn,1);
    % mat(keptvox(sig))=peaks(sig);
    % nii=mat2nii(mat);
    % save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_fdr_peaks.nii' ]);
    %
    % mat=nan(voxn,1);
    % mat(keptvox)=peaks;
    % nii=mat2nii(mat);
    % save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peaks.nii' ]);
    %
    %
    mat=nan(voxn,1);
    mat(keptvox)=peakLags;
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReverse_peakLags.nii' ])
    
    % mat=nan(voxn,1);
    % mat(keptvox(sig))=peakLags(sig);
    % nii=mat2nii(mat);
    % save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReverse_fdr_peakLags.nii' ])
    %
    %  peakLags_bined=discretize(peakLags,[-15 -10 -6 -3 0 1 4  7  11 16]);
    % mat=nan(voxn,1);
    % mat(keptvox)=peakLags_bined;
    % nii=mat2nii(mat);
    % save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peakLagsBined.nii' ]);
    
    peakLags_bined=discretize(peakLags,[-15 -1 0 1 2 16]);
    mat=nan(voxn,1);
    mat(keptvox)=peakLags_bined;
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_peakLags5bins.nii' ]);
    
    mat(keptvox(sig==0))=NaN;;
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_pfdr_peakLags5bins.nii' ]);
    
    peakLags_bined=discretize(peakLags,[-15 -1 0 1 2 16]);
    mat=nan(voxn,1);
    mat(keptvox(abs(npeaks)<peaks | isnan(npeaks)))=peakLags_bined(abs(npeaks)<peaks | isnan(npeaks));
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReversed_maxpos_peakLags5bins.nii' ]);
    
end
