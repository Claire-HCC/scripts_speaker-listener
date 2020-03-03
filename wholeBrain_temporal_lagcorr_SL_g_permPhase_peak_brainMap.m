% loc='cluster';
set_parameters;
timeUnit='tr' ;
lags_tested={-10:10,  -30:30};
permN=1000;

for ei=1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load(sprintf('%s/%s/fmri/temporal/lagcorr/%s/wholeBrain/SL_g/perm/lag%d-%d_permPhase',expdir,exp,timeUnit,min(lags),max(lags)),'r','keptvox');
        r_perm=r;
        
   %     load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'sig_fdr','keptvox');
    %    regression_mask=sig_fdr;
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r','lags','keptT','keptvox');
        r=r;
        
        p=mean(r_perm>r,3);
        %    p(regression_mask==0,:)=NaN;
        
        [peak peak_lagi]=max(r,[],2);
        peakLags=(lags(peak_lagi))';
        %  peakLags(regression_mask==0)=NaN;
        
        pfdr=nan(size(p(:)));
        if sum(~isnan(p))>0
            [~,~,pfdr(~isnan(p))]=fdr(p(~isnan(p)));
            pfdr=reshape(pfdr,size(p));
        end
        
        r_temp=r;
        r_temp(pfdr>.05)=NaN;
        peak_pfdr=nan(size(peak));
        peakLags_pfdr=peak_pfdr;
        vis=find(sum(~isnan(r_temp),2)~=0);
        [peak_pfdr(vis) peak_lagi]=max(r_temp(vis,:),[],2);
        peakLags_pfdr(vis)=lags(peak_lagi);
        
        mat=nan(voxn,1);
        mat(keptvox)=peakLags+0.00000001;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=peakLags_pfdr;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peakLags_pfdr.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=peak;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak.nii']);
        
        mat=nan(voxn,1);
        mat(keptvox)=peak_pfdr;
        nii=mat2nii(mat);
        save_nii(nii,[expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags))  '_peak_pfdr.nii']);
        
        save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/wholeBrain/SL_g/lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'keptT','lags','keptvox','p','pfdr','peak','peakLags','peak_pfdr','peakLags_pfdr','-v7.3');
        clear peak peak_lagi peakLags
    end
end

