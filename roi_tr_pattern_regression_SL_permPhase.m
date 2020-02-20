function roi_tr_pattern_regression_SL_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags_tested={-10:10,  -30:30};
% set random seed. Otherwise, phase_rand2 will always get the same result.
permN=1000;

rname=rnames{ri};
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        b=nan([1 length(lags)+1 permN]);
        r2=nan([1 permN ]);
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat'],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/speaker_' rname '.mat'],'data');
            data_perm=zscore(data,0,2);
            
            [roi_voxn, tn,~]=size(gdata);
            
            keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            g=nanmean(gdata,3);
            y=g;
            y=y(:,keptT);
            y=y(:);
            
            for perm=1:permN;
                for li=1:length(lags);
                    X(:,:,li)=data_perm(:,keptT+lags(li),perm);
                end
                
                X=reshape(X,roi_voxn*length(keptT),length(lags));
                % centralized X
                X=X-mean(X);
                
                % add an constant
                [b(1,:,perm),~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
                r2(1,perm)=stats(1);
                
                clear X
            end
        end
        
        save(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/SLg/perm/regression_SL_lag%d-%d_permPhase_%s',expdir,exp,timeUnit,froidir,min(lags),max(lags),rname),'b','r2','lags','keptT');
        clear b r2 keptT rp
    end
end

