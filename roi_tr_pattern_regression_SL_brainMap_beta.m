set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

for ei=[1 ];%[1:4];;%
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats.mat'],'b_real','b_null''b_sig_fdr','b_sig_fwe','r2_sig_fdr','r2_sig_fwe','lags');
        
        b=b_real;
       b(r2_sig_fwe==0,:)=NaN;
        [bPeak, bPeakLagi]=max(b');
        bPeakLags(:,1)=(lags(bPeakLagi));
        
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(~isnan(bPeakLags)),'UniformOutput',0);
        roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
        nii=roiTable2wholeBrainNii_mor([roi_ids, bPeakLags(~isnan(bPeakLags))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_bPeakLags.nii']);
    end
    
end

        