
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '\roi_mask\mor\' 'roi_id_region.mat'],'roi_table');
win_width=100;
win_step=50;

for ei=1%:2;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag-5-5' ],'couplingz','rnames','keptT');
    
    
    cp_sl=couplingz;
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag-5-5' ],'couplingz','rnames','keptT');
    cp_ll=nanmean(couplingz,3);
    
    for ri=1:length(rnames);
       
        winN=floor((length(keptT)-win_width+1)/win_step)+1;
        for wi=1:winN;
            win_s=(wi-1)*win_step+1;
            win_e=win_s+win_width-1;
            [herd(ri,wi) p(ri,wi)]=corr(cp_sl(ri,keptT(win_s:win_e))',cp_ll(ri,keptT(win_s:win_e))');
        end
    end
    herd([3 4 35 54 ],:)
    sig=p([3 4 35 54 ],:)
    sig=fdr0(sig(:),0.05);
     sig=reshape(sig,length(sig)/winN,winN)
  %  sig=reshape(sig,size(p));
  %  rnames(sum(sig,2)>0)
end


