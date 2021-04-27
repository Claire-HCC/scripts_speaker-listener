function roi_temporal_regression_SL_g_permPhase(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags_tested={-10:-4. 4:10,0 ,-10:10, -40:40};
% set random seed. Otherwise, phase_rand2 will always get the same result.
permN=1000;
rname=rnames{ri};

for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        b=nan([1 length(lags)+1 permN ]);
        r2=nan([1 permN ]);
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat'],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            gdata(:,:,subjects_excluded{ei})=NaN;
            data=nanmean(data,1);
            gdata=nanmean(gdata,1);
            
            [~,tn,subjn]=size(gdata);
            roi_voxn=1;
            
            keptT_s= max(-min(lags),0)+crop_start+1;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            y=nanmean(gdata(:,keptT,:),3);
            y=y(:);
            
            for perm=1:permN;
                
                % average all voxels within the roi
                data_perm=nan([1 tn]);
                data_perm((crop_start+1):tn)=(phase_rand2(data(:,(crop_start+1):tn)',1))';
                
                for li=1:length(lags);
                    X(:,:,li)=data_perm(:,keptT+lags(li));
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
        
        save(sprintf('%s/%s/fmri/temporal/regression/%s/roi/%s/SL_g/perm/lag%d-%d_permPhase_%s',expdir,exp,timeUnit,froidir,min(lags),max(lags),rname),'b','r2','keptT');
        clear b r2 keptT rp
    end
end
