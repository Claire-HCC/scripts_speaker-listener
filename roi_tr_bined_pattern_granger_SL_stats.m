clear all
close all
% loc='cluster';
set_parameters;
exp='merlin';
timeUnit='tr' ;
froidir='mor';
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=0;
binSize_tested=[10 30]; % tr;
binStep=1;
lags_tested={-10:-1,-10:-4};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

for binSizei=1:length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'p','b_l_perm','b_ls_perm','b_s_perm','F_l_perm','F_ls_perm','F_s2l_perm','F_l2l_perm','F_s_perm','r2_ls_perm','r2_l_perm','r2_s_perm','keptT','rnames');
        keptT_perm=keptT;
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))  ],'F_s2l','keptT','rnames');
        keptT_real=keptT;
        keptT=min(intersect(keptT_perm,keptT_real)):min(size(F_s2l_perm,2),size(F_s2l,2));
        
        p=nan(size(F_s2l));
        p(:,keptT)=(sum(F_s2l_perm(:,keptT,:)>F_s2l(:,keptT),3))/size(F_s2l,1);
        figure;; imagesc(p==0);
        
        keptT=keptT_perm;
          save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'p','b_l_perm','b_ls_perm','b_s_perm','F_l_perm','F_ls_perm','F_s2l_perm','F_l2l_perm','F_s_perm','r2_ls_perm','r2_l_perm','r2_s_perm','keptT','rnames');

    end
end
