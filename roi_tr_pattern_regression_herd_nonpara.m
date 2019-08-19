
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '\roi_mask\mor\' 'roi_id_region.mat'],'roi_table');

for ei=2%:2;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag-5-5' ],'couplingz','rnames','keptT');
    cp_sl=couplingz;
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag-5-5_perm' ],'couplingz_perm','rnames','keptT');
    cp_sl_perm=couplingz_perm;
    iters=size(couplingz_perm,3);
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag-5-5' ],'couplingz','rnames','keptT');
    cp_ll=nanmean(couplingz,3);
    
    for ri=1:length(rnames);
        [herd(ri,1) p(ri,1)]=corr(cp_sl(ri,keptT)',cp_ll(ri,keptT)');
        
        for iter=1:iters;
            [herd_perm(ri,1,iter) ]=corr(cp_sl_perm(ri,keptT,iter)',cp_ll(ri,keptT)');
            
        end
    end
    
    sig=fdr0(p(:),0.05);
    rnames(sum(sig,2)>0)
end

