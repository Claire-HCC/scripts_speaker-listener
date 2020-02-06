% bined r2 vs std between subjects
clear all;
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

% -4 for merlin, -3 for sherlock
lags=-6:-3;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

for ei=3;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/vox_regression_LLminus1L_lag0-0'],'r','gdata');
    % cp_ll=squeeze(nanmean(r,4));
    cp_ll=repmat(std(gdata,[],3)',1,18);
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/vox_regression_SLminus1L_lag' num2str(min(lags)) '-' num2str(max(lags))],'r','voxind','keptvox','keptT');
    cp_sl=squeeze(abs(r));
    
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/vox_regression_SLfake_lag' num2str(min(lags)) '-' num2str(max(lags))],'r','gdata');
    cp_slf=squeeze(abs(r));
    

    for sp=1:size(r,3);
        herd(sp,1)=corr_col(cp_sl(keptT,sp),cp_ll(keptT,sp));
        herd_null(sp,1)=corr_col(cp_slf(keptT,sp),cp_ll(keptT,sp));
    end
    
    
    herdz=real(atanh(herd));
    herdz_m=mean(herdz,2);
    herdz_null=real(atanh(herd_null));
    for vi=1:size(r,1);
        [~,herd_p(vi,1)]=ttest(herdz-herdz_null,0,'Tail','right');
    end
    
    % test herd effect within regions showing significant coupling
    
    herd_sig_fdr=fdr0(herd_p,0.05);
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/herd_vox_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'lags','herd','herd_p','herd_sig_fdr','herd_null','Y','y');
    
    % clear herd herd_null herd_p herd_sig_fdr
end
