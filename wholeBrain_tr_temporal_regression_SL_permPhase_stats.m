clear all
close all

% loc='cluster';
set_parameters;
timeUnit='tr' ;
lags=-10:10;

for ei=1:4;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain//SLg/perm/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase.mat' ],'r2_perm','b_perm','keptvox');
    load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain//SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat' ],'r2','b','keptvox');

    p=sum(r2_perm>r2(:,1),2)/size(r2_perm,2);
    p_fdr=nan(size(p));
    ris=find(~isnan(r2(:,1)));
    [~,~,p_fdr(ris)]=fdr(p(ris));
    sig_fdr=p_fdr<.05;
    
    p_fwemax=sum(repmat(max(r2_perm),length(r2),1)>r2,2)/(size(r2_perm,2));
    p_fwemax(isnan(r2))=NaN;
    
    mat=nan(voxn,1);
    mat(keptvox( p_fwemax<.05))=r2(p_fwemax<.05);
    nii=mat2nii(mat);
    save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_stats.mat' ],'r2','sig_fdr','p','p_fdr','p_fwemax','r2_perm','b_perm','b','-v7.3');
    
    mat=nan(voxn,1);
    mat(keptvox( p_fwemax<.05))=r2(p_fwemax<.05);
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_r2_fwemax.mat' ]);
    
    mat=nan(voxn,1);
    mat(keptvox( p_fdr<.05))=r2(p_fdr<.05);
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/temporal_regression/' timeUnit '/wholeBrain/SLg/lag' num2str(min(lags)) '-' num2str(max(lags)) '_permPhase_r2_fdr.mat' ])

end


