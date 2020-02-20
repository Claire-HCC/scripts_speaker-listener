
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};
permN=1000;

for ei=1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'sig_fdr');
        regression_mask=sig_fdr;
        
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        p=sum(r_perm>r,3)/permN;
        p(regression_mask==0,:)=NaN;
        
        
        p_fdr=nan(size(p(:)));
        if sum(~isnan(p))>0
        [~,~,p_fdr(~isnan(p))]=fdr(p(~isnan(p)));
        p_fdr=reshape(p_fdr,size(p));
        end
        
        onsets=nan(length(rnames),1);
        for ri=1:length(rnames);
            if sum(p_fdr(ri,:)<.05)>0;
                [lagi]=min(find(p_fdr(ri,:)<.05));
                onsets(ri)=lags(lagi);
            end
        end
        table(onsets(~isnan(onsets)),rnames(~isnan(onsets)))
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( onsets)),   onsets(~isnan( onsets))]);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_betaPeak_neg.nii']);
            
        save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','p_fdr','onsets');
    end
end

