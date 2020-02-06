clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

lags=-10:10;

% % SL
for ei=1:4;%1:4;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F_sl','r2_sl','p_sl','keptvox','b_sl');
    
    
    couple_sl_sig_fwe=(p_sl<(0.05/length(F_sl)));
    couple_sl_sig_fdr=fdr0(p_sl,0.05);
    
    %     save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_sl_sig_fdr','couple_sl_sig_fwe','p_sl','lags','keptvox','r2_sl','F_sl');
    %
    %     mat=zeros(voxn,1);
    %     mat(keptvox(couple_sl_sig_fdr==1))=F_sl(couple_sl_sig_fdr==1);
    %     nii=mat2nii(mat);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_F_fdr' ]);
    %
    %     mat=zeros(voxn,1);
    %     mat(keptvox(couple_sl_sig_fwe==1))=F_sl(couple_sl_sig_fwe==1);
    %     nii=mat2nii(mat);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_F_fwe' ]);
    %
    %     mat=zeros(voxn,1);
    %     mat(keptvox(couple_sl_sig_fdr==1))=r2_sl(couple_sl_sig_fdr==1);
    %     nii=mat2nii(mat);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_r2_fdr' ]);
    %
    %     mat=zeros(voxn,1);
    %     mat(keptvox(couple_sl_sig_fwe==1))=r2_sl(couple_sl_sig_fwe==1);
    %     nii=mat2nii(mat);
    %     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_r2_fwe' ]);
    
    % find peak lags
    b_temp=b_sl(:,2:end);
    [v col]=max(b_temp(:,2:end),[],2);
    bpt=nan(size(col));
    bpt(:,1)=lags(col);
    mat=zeros(voxn,1);
    mat(keptvox)=bpt;
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/coulpingPeakLag_sl_'  num2str(min(lags)) '-' num2str(max(lags)) '.nii' ]);
    
    
        % focus on voxels that showed significant SL coupling
    b_temp=b_sl(:,2:end);
    b_temp(couple_sl_sig_fdr~=1)=0;
    [v col]=max(b_temp(:,2:end),[],2);
    col(sum(couple_sl_sig_fdr,2)==0)=nan;
    bpt=nan(size(col));
    bpt(~isnan(col),1)=lags(col(~isnan(col)));
    mat=zeros(voxn,1);
    mat(keptvox(couple_sl_sig_fdr==1))=bpt(couple_sl_sig_fdr==1);
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/coulpingPeakLag_sl_'  num2str(min(lags)) '-' num2str(max(lags)) '_fdr.nii' ]);
    
    
end

%% LL
% for ei=1:4;
%     exp=experiments{ei};
%     load([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2_ll','p_ll','keptvox');
%
%     couple_ll_sig_fwe=(p_ll<(0.05/size(r2_ll,1)));
%     couple_ll_sig_fwe=(squeeze(sum(couple_ll_sig_fwe==0,3))==0);
%
%     for s=1:size(r2_ll,3);
%         couple_ll_sig_fdr(:,s)=fdr0(squeeze(p_ll(:,1,s)),0.05);
%     end
%     couple_ll_sig_fdr=(sum(couple_ll_sig_fdr==0,2)==0);
%
%     r2_ll_g=mean(r2_ll,3);
%     save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_LL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_stats' ],'couple_ll_sig_fdr','couple_ll_sig_fwe','lags','keptvox','r2_ll_g');
%
%     mat=zeros(voxn,1);
%     mat(keptvox(sum(couple_ll_sig_fdr==0,2)==0))=r2_ll_g(couple_ll_sig_fdr==1);
%     nii=mat2nii(mat);
%     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_LL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_r2_fdr' ]);
%
%     mat=zeros(voxn,1);
%     mat(keptvox(couple_ll_sig_fwe==1))=r2_ll_g(couple_ll_sig_fwe==1);
%     nii=mat2nii(mat);
%     save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_LL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_r2_fwe' ]);
%
%     clear couple_ll_sig_fdr
% end
