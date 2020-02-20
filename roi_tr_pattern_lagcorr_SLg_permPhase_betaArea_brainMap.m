
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};
permN=1000;

for ei=3;%2:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        p=sum(r_perm>r,3)/permN;
        
        % listner precede area
        r_temp=r;
        r_temp(p>.05)=NaN;
        r_temp(:,lags<0)=NaN;
        r_temp_norm=r_temp./repmat(nansum(r_temp,2),1,length(lags));
        a=lags.*r_temp_norm;
        lag_pos(:,1)=nansum(a(:,lags>0),2);
        lag_pos(sum(~isnan(r_temp_norm),2)==0)=NaN;
        
        % speaker precedes area
        r_temp=r;
        r_temp(p>.05)=NaN;
        r_temp(:,lags>0)=NaN;
        r_temp_norm=r_temp./repmat(nansum(r_temp,2),1,length(lags));
        a=lags.*r_temp_norm;
        lag_neg(:,1)=nansum(a(:,lags<0),2);
        lag_neg(sum(~isnan(r_temp_norm),2)==0)=NaN;
        
        % both areas
        r_temp=r;
        r_temp(p>.05)=NaN;
        r_temp_norm=r_temp./repmat(nansum(r_temp,2),1,length(lags));
        a=lags.*r_temp_norm;
        lag_both(:,1)=nansum(a,2);
        lag_both(sum(~isnan(r_temp_norm),2)==0)=NaN;
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( lag_neg)),   lag_neg(~isnan(  lag_neg))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_betaArea_neg.nii']);
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( lag_pos)),   lag_pos(~isnan(  lag_pos))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_betaArea_pos.nii']);
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( lag_both)),   lag_both(~isnan(  lag_both))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_betaArea_both.nii']);
        
        save([expdir '/' exp '/fmri/pattern_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'rnames','r','lags','keptT','r_perm','p','lag_pos','lag_neg','lag_both');
    end
end

