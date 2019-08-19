
clear all;
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';
rnames=dir([expdir '/roi_mask/'  froidir '/mat/*.mat']);
rnames=strrep({rnames.name},'.mat','');

tic % 15 min
win_width=25;
win_step=1;
lags=-20:20;

for ei=1%:2;
    exp=experiments{ei};
    
    for ri=44;%1:length(rnames);
        rname=rnames{ri};
        fl=[expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname  ];
        
        if exist([fl '.mat'])>0;
            load([expdir '/' exp '/fmri/isc/' timeUnit '/roi/' froidir '/iscz_SL_' rname ],'isc','lags');
            isc_sl=isc;
            load([expdir '/' exp '/fmri/isc/' timeUnit '/roi/' froidir '/iscz_LL_' rname ],'isc','lags');
            isc_ll=isc;
            
             load([expdir '/' exp '/fmri/isc/' timeUnit '/roi/' froidir '/iscz_SL_' rname  '_perm'],'isc_perm','lags');
             isc_sl_perm=isc_perm;
             load([expdir '/' exp '/fmri/isc/' timeUnit '/roi/' froidir '/iscz_LL_' rname '_perm'],'isc_perm','lags');
             isc_ll_perm=isc_perm;
             
            [isc_sl_bestT_r,isc_sl_bestT]=max(nanmean(nanmean(isc_sl,3),2));
            [isc_ll_bestT_r,isc_ll_bestT]=max(nanmean(nanmean(isc_ll,3),2));
           
             [isc_sl_perm_bestT_r,isc_sl_perm_bestT]=max(nanmean(nanmean(isc_sl_perm,3),2));
             sum(isc_sl_bestT_r<isc_sl_perm_bestT_r)/100;
             
            isc_sl_m=nanmean(isc_sl,3);
             isc_sl_perm_m=nanmean(isc_sl_perm,3);
            isc_ll_m=nanmean(isc_ll,3);
            
            herd=corr(isc_sl_m(isc_sl_bestT,:)',isc_ll_m(lags==0,:)');
            
            iters = size(isc_sl_perm,4);
            for i=1:iters;
                isc_sl_perm_m_bestT(:,i)=isc_sl_perm_m(isc_sl_perm_bestT(i),:,:,i);
            end
            
            herd_perm=corr_col(isc_sl_perm_m_bestT,repmat(isc_ll_m(lags==0,:)',1,iters));
            
            isc=0.5*log((1+isc_SL)./(1-isc_SL));
            save([expdir '/' exp '/fmri/isc/' timeUnit '/roi/' froidir '/iscz_SL_' rname ],'isc','lags');
            isc=0.5*log((1+isc_LL)./(1-isc_LL));
            save([expdir '/' exp '/fmri/isc/' timeUnit '/roi/' froidir '/iscz_LL_' rname  ],'isc','lags');
            
        end
    end
end


toc
