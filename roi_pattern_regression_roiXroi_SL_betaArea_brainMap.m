
% loc='cluster';
clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, -10:-4, -20:-4, -30:-4, -10:-1};

for ei=1:4;%
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        load(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLeach/regression_SLeach_roiXroi_lag%d-%d',expdir,exp,timeUnit,froidir,min(lags),max(lags)),'b','r2','r','p','F');
        
        blag_pos=nan([length(rnames) length(rnames)]);
        blag_neg=blag_pos;
        
        for ri_s=1:length(rnames);
            for ri_l=1:length(rnames);
                
                % within roi showing signifcant difference in sl and ll beta
                % profiles.
                %    if sig_betaClass(ri)==1;
                
                b_s_temp=squeeze(b(ri_s,ri_l,2:end,:))';
                
                % find the biggest peak within lag<=0
                [~,p,~,stats]=ttest(b_s_temp,0,'tail','right');
                
                b_norm=(nanmean(b_s_temp,1)/sum(nanmean(b_s_temp(:,lags<0 & p<.05),1)));
                a=lags.*b_norm.*(p<.05);
                blag_neg(ri_s,ri_l)=nansum(a(lags<0));
                
                b_norm=(nanmean(b_s_temp,1)/sum(nanmean(b_s_temp(:,lags>0 & p<.05),1)));
                a=lags.*b_norm.*(p<.05);
                blag_pos(ri_s,ri_l)=nansum(a(lags>0));
            end
        end
        lag_neg(lag_neg==0)=NaN;
        lag_pos(lag_pos==0)=NaN;
           lag_both(lag_both==0)=NaN;
        
        save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_roiXroi_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bLagArea'  ],'rnames','blag_pos','blag_neg');
        
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan(lag_neg)),  lag_neg(~isnan( lag_neg))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_roiXroi_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_blag_neg.nii']);
        nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( blag_pos)),   blag_pos(~isnan(  blag_pos))]);
        save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_roiXroi_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_blag_pos.nii']);
    end
end









