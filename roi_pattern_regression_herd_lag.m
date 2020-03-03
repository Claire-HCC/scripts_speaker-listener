
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '\roi_mask\mor\' 'roi_id_region.mat'],'roi_table');

lags=-20:20;
for ei=1%:2;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag-5--1' ],'couplingz','rnames','keptT');
    cp_sl=couplingz;

    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag0-0' ],'couplingz','rnames','keptT');
    cp_ll=nanmean(couplingz,3);
    
    for ri=1:length(rnames);
        [herd(ri,:) ]=lagcorr(cp_ll(ri,keptT)',cp_sl(ri,keptT)',lags);
    end
    
    sig=fdr0(p(:),0.05);
    sig=reshape(sig,size(p));
    rnames(sum(sig,2)>0)
end

