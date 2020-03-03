
set_parameters; 
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -30:30};
permN=1000;

for ei=1:4;%2:4;
    exp=experiments{ei};
    
    for lagi=1;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/perm/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase' ],'r');
        r_perm=r;
        
        load([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');
        p=sum(r_perm>r,3)/permN;
        
        % listner precede area
        r_temp=r;
        r_temp(p>(.05/(length(rnames)*length(lags))))=NaN;
        
        lag_pos=nan([length(rnames) 1]);
        ris=find(nansum(r_temp(:,lags>0),2)>0);
        if ~isempty(ris);
            [~,lagi]=max(r_temp(:,lags>0),[],2);
            lag_pos(ris)=lags(lagi(ris)+sum(lags<=0));
        end
        
        lag_neg=nan([length(rnames) 1]);
        ris=find(nansum(r_temp(:,lags<0),2)>0);
        if ~isempty(ris);
            [~,lagi]=max(r_temp(:,lags<0),[],2);
            lag_neg(ris)=lags(lagi(ris));
        end
        
        lag_both=nan([length(rnames) 1]);
        ris=find(nansum(r_temp,2)>0);
        if ~isempty(ris);
            [~,lagi]=max(r_temp,[],2);
            lag_both(ris)=lags(lagi(ris));
        end
        
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( lag_neg)),   lag_neg(~isnan(  lag_neg))]);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_betaPeak_neg.nii']);
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( lag_pos)),   lag_pos(~isnan(  lag_pos))]);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_betaPeak_pos.nii']);
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( lag_both)),   lag_both(~isnan(  lag_both))]);
        save_nii(nii,[expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_betaPeak_both.nii']);
        
        save([expdir '/' exp '/fmri/temporal_lagcorr/' timeUnit '/roi/' froidir '/SLg/SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_betaPeaks' ],'rnames','r','lags','keptT','r_perm','p','lag_pos','lag_neg','lag_both');
    end
end

