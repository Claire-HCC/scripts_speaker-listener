
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

tic % 15 min
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

for ei=1:4;
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for lagi=6;%1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
        b_real=b(:,2:end);
        F_real=F;
        r2_real=r2;
        r2_byTime_real=r2_byTime;
        tn=size(r2_byTime,2);
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rnames{1} '.mat' ],'gdata');
        listenersN=size(gdata,3);
        
        for s=1:listenersN;
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/perm/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' num2str(s)],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
            b_null(:,:,s)=b(:,2:end);
            F_null(:,s)=F;
            r2_byTime_null(:,:,s)=r2_byTime;
            r2_null(:,s)=r2;
        end
        
        r2_p=(sum(r2_null>r2,2)/listenersN);
        r2_sig_fdr=(fdr0(r2_p,0.05)==1);
        r2_sig_fwe=nan(size(r2_p));
        r2_sig_fwe=(r2_p<(0.05/length(rnames)));
        
        r2_byTime_p=(sum(r2_byTime_null>r2_byTime_real,3)/listenersN);
        
        r2_byTime_sig_fdr=nan(size(r2_byTime_p));
        for ri=find(r2_sig_fdr==1);
            r2_byTime_sig_fdr(ri,keptT)=(fdr0(r2_byTime_p(ri,keptT),0.05)==1);
        end
        
        b_p=(sum(b_null>b_real,3)/listenersN);
        b_sig_fdr=nan(size(b_p));
        b_sig_fwe=nan(size(b_p));
        for ri=find(r2_sig_fdr==1);
            b_sig_fdr(ri,:)=(fdr0(b_p(ri,:),0.05)==1);
            b_sig_fwe(ri,:)=((b_p(ri,:)<(0.05/length(lags)))==1);
        end
        
        save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats.mat'],'b_real','b_null','b_p','b_sig_fdr','b_sig_fwe','r2_p',...
            'r2_sig_fdr','r2_sig_fwe','r2_byTime_p','r2_byTime_sig_fdr','lags','rnames','keptT','r2_byTime_null','r2_byTime_real');
        
        
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(r2_sig_fwe==1 ),'UniformOutput',0);
        roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
        nii=roiTable2wholeBrainNii_mor([roi_ids, r2(r2_sig_fwe==1)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_r2_tmap_fwe.nii']);
        
        roi_table_inds=cellfun(@(x) strmatch(x,roi_table.region,'exact'),rnames(r2_sig_fdr==1 ),'UniformOutput',0);
        roi_ids=cell2mat(roi_table.id(cell2mat(roi_table_inds)));
        nii=roiTable2wholeBrainNii_mor([roi_ids, r2(r2_sig_fdr==1)]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_r2_tmap_fdr.nii']);
        
        clear b_real b_null b_p b_t b_sig_fdr b_sig_fwe r2_p r2_t r2_sig_fdr r2_sig_fwe r2_byTime_p r2_byTime_t r2_byTime_sig_fdr r2_byTime_null r2_byTime_real F_null r2_null
    end
end



