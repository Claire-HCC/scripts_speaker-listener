% loc='cluster';
set_parameters;
timeUnit='tr' ;
lags_tested={-10:10,  -30:30};
permN=1000;
load([expdir '/roi_mask/gray_matter_mask.mat'],'roimask');
mask_gray=find(roimask==1);

for ei=[3:4];%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load(sprintf('%s/%s/fmri/temporal_lagcorr/%s/wholeBrain/SLg/perm/lag%d-%d_permPhase',expdir,exp,timeUnit,min(lags),max(lags)),'r','keptvox');
        r_perm=r(ismember(keptvox,mask_gray),:,:);
        
        load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'sig_fdr','keptvox');
        regression_mask=sig_fdr(ismember(keptvox,mask_gray));
        
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptT','keptvox');
        r=r(ismember(keptvox,mask_gray),:);
        
        keptvox=keptvox(ismember(keptvox,mask_gray));
        p=mean(r_perm>r,3);
    %    p(regression_mask==0,:)=NaN;
        
        [peak peak_lagi]=max(r,[],2);
        peakLags=(lags(peak_lagi))';
      %  peakLags(regression_mask==0)=NaN;
        
        p_fdr=nan(size(p(:)));
        if sum(~isnan(p))>0
            [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
            p_fdr=reshape(p_fdr,size(p));
        end
        
        onsets_pfdr=nan(length(keptvox),1);
        onsets_p=nan(length(keptvox),1);
        for vi=1:length(keptvox);
            if sum(p_fdr(vi,:)<.05)>0;
                [lagi]=(find(p_fdr(vi,:)<.05));
                onsets_pfdr(vi)=lags(min(lagi));
            end
            
            if sum(p(vi,:)<.05)>0;
                [lagi]=(find(p(vi,:)<.05));
                onsets_p(vi)=lags(min(lagi));
            end
        end
        %        table(onsets(~isnan(onsets)),keptT(~isnan(onsets)))
        
        [~,~,onsetsOrder]=unique(onsets_pfdr);
        onsetsOrder(isnan(onsets_pfdr))=NaN;
        
        mat=nan(voxn,1);
        mat(keptvox)=peakLags+0.00000001;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags.nii']);
        
        mat=nan(voxn,1);
        peakLags_pfdr=peakLags+0.00000001;
        peakLags_pfdr(isnan(onsets_pfdr))=NaN;
        mat(keptvox)=peakLags_pfdr;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags_pfdr.nii']);
        
        mat(keptvox)=peakLags+0.00000001;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=peak;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak.nii']);
        
        peak_pfdr=peak;
        peak_pfdr(isnan(onsets_pfdr))=NaN;
        mat=nan(voxn,1);
        mat(keptvox)=peak_pfdr;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak_pfdr.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=onsets_pfdr+0.00000001;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsets_pfdr.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=onsets_p+0.00000001;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsets_p.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=onsetsOrder;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsetsOrder_pfdr.nii']);
        
        save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'keptT','lags','keptvox','p','p_fdr','onsets_p','onsets_pfdr','peak','peakLags','-v7.3');
        clear peak peak_lagi peakLags
    end
end

