
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};

for ei=1:4;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'sig_fdr');
        regression_mask=sig_fdr;
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        p=mean(r_perm>r,3);
      %  p(regression_mask==0,:)=NaN;
        p(isnan(r(:,1)))=NaN;
        
        [peak peakLagi]=max(r,[],2);
        peakLags=(lags(peakLagi))';
       % peakLags(regression_mask==0,:)=NaN;
        peakLags(isnan(r(:,1)),:)=NaN;
        
        p_fdr=nan(size(p(:)));
        if sum(~isnan(p))>0
            [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
            p_fdr=reshape(p_fdr,size(p));
        end
        
        onsets_pfdr=nan(length(rnames),1);
        onsets1=nan(length(rnames),1);
        onsets2=nan(length(rnames),1);
        for ri=1:length(rnames);
            if sum(p_fdr(ri,:)<.05)>0;
                [lagi]=(find(p_fdr(ri,:)<.05));
                onsets_pfdr(ri)=lags(min(lagi));
                
                try
                    onsets1(ri)=min(intersect(lags(lagi),lags(lags<0)));end
                try
                    onsets2(ri)=min(intersect(lags(lagi),lags(lags>=0)));end
            end
        end
        table(onsets_pfdr(~isnan(onsets_pfdr)),rnames(~isnan(onsets_pfdr)))
        
        [~,~,onsetsOrder]=unique(onsets_pfdr);
        onsetsOrder(isnan(onsets_pfdr))=NaN;
        
        %   delete([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/*nii']);
        % without adding 0.0001, zeros lag is not displayed with xjview
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( onsets_pfdr)),   onsets_pfdr(~isnan( onsets_pfdr))+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsets_pfdr.nii']);
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( onsetsOrder)),   onsetsOrder(~isnan( onsetsOrder))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsetsOrder_pfdr.nii']);
        
        peak_pfdr=peak;
        peak_pfdr(isnan(onsets_pfdr))=NaN;
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peak_pfdr)),   peak_pfdr(~isnan( peak_pfdr))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/'  froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak_pfdr.nii']);
        
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags)),   peakLags(~isnan(peakLags))+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/'  froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags.nii']);
        
        peakLags_pfdr=peakLags;
        peakLags_pfdr(isnan(onsets_pfdr))=NaN;
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( peakLags_pfdr)),   peakLags_pfdr(~isnan( peakLags_pfdr))+0.00000001]);
        nii.img(1,1,1)=1;
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/'  froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags_pfdr.nii']);
        
        save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','p_fdr','onsets_pfdr','onsetsOrder','peak','peakLags','peakLags_pfdr');
    end
end

