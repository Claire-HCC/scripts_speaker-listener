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

for ei=[1];;%
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'r2','sig_fdr','p','p_fdr','p_fwemax','r2_perm','b_perm','b');
        
        blag_pos=nan([length(rnames) 1]);
        blag_neg=blag_pos;
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            
            % within roi showing signifcant difference in sl and ll beta
            % profiles.
            %    if sig_betaClass(ri)==1;
            b_temp=b(ri,2:end);
            b_perm_temp=squeeze(b_perm(ri,2:end,:))';
            
            % find the biggest peak within lag<=0
            p=sum(b_perm_temp>b_temp,1)/size(b_perm,3);
            b_norm=b_temp/sum(b_temp(lags<0 & p<.05));
            a=lags.*b_norm.*(p<.05);
            blag_neg(ri,1)=nansum(a(lags<0));
            
            b_norm=b_temp/sum(b_temp(lags<0 & p>.05));
            a=lags.*b_norm.*(p<.05);
            blag_pos(ri,1)=nansum(a(lags>0));
        end
    end
    blag_neg(blag_neg==0)=NaN;
    blag_pos(blag_pos==0)=NaN;
    save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bLagArea'  ],'rnames','blag_pos','blag_neg');
    
    nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( blag_neg)),   blag_neg(~isnan(  blag_neg))]);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_blag_neg.nii']);
    nii=roiTable2wholeBrainNii_mor([roi_ids(~isnan( blag_pos)),   blag_pos(~isnan(  blag_pos))]);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))  '_blag_pos.nii']);
end


