
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};
permN=1000;
binSize=40;

for ei=3;%1:4;%[1 2 4];%2:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'sig_fdr');
        regression_mask=sig_fdr;
        
         load([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/perm/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags');
        p=sum(r_perm>r,4)/permN;
        p(regression_mask==0,:,:)=NaN;
        p(isnan(r))=NaN;
        [~,tn,~]=size(r);
        
        p_fdr=nan(size(p(:)));
        if sum(~isnan(p(:)))>0
            [~,~,p_fdr(~isnan(p(:)))]=fdr(p(~isnan(p(:))));
            p_fdr=reshape(p_fdr,size(p));
        end
        
        onsets=nan(length(rnames),tn);
        for ri=1:length(rnames);
            for t=1:tn;
                if nansum(squeeze(p_fdr(ri,t,:)<.05))>0;
                    [lagi]=min(find(p_fdr(ri,t,:)<.05));
                    onsets(ri,t)=lags(lagi);
                end
            end
        end
        
        [~,~,onsetsOrder]=unique(onsets(:));
        onsetsOrder(isnan(onsets(:)))=NaN;
        onsetsOrder=reshape(onsetsOrder,size(onsets));
        
        %         nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( onsets)),   onsets(~isnan( onsets))]);
        %         save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsets.nii']);
        %
        %         nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( onsetsOrder)),   onsetsOrder(~isnan( onsetsOrder))]);
        %         save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_onsetsOrder.nii']);
        %
        save([expdir '/' exp '/fmri/pattern_lagcorr_bined/' timeUnit '/roi/' froidir '/SLg/SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_onsets' ],'rnames','r','lags','keptT','r_perm','p','p_fdr','onsets','onsetsOrder');
    end
end

