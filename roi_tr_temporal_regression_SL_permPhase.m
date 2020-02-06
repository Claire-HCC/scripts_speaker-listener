function roi_tr_temporal_regression_SL_permPhase(perm)

% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags=-10:10;
% set random seed. Otherwise, phase_rand2 will always get the same result.
rng(perm)

for ei=1:4;
    exp=experiments{ei};
    
    b=nan([length(rnames) length(lags)+1 ]);
    r2=nan([length(rnames) 1 ]);
    
    for ri=1:size(rnames);
        rname=rnames{ri};
        
        disp(ri)
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat'],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat'],'data');
            
            % average all voxels within the roi
            gdata=nanmean(gdata,1);
            data=nanmean(data,1);
            
            roi_voxn=size(gdata,1);
            tn=size(gdata,2);
            
            keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            data_perm=nan(size(data))';
            if ri==1;
                [data_perm((crop_start+1):end,:), rp(:,perm)]=phase_rand2(data(:,(crop_start+1):end)',1);
            else;
                [data_perm((crop_start+1):end,:)]=phase_rand2(data(:,(crop_start+1):end)',1,rp(:,perm));
            end
            data_perm=data_perm';
            
            g=nanmean(gdata,3);
            y=g;
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data_perm(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            % centralized X
            X=X-mean(X);
            
            % add an constant
            [b(ri,:),~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
            r2(ri,1)=stats(1);
            
            clear X
        end
    end
    
    save(sprintf('%s/%s/fmri/temporal_regression/%s/roi/%s/SLg/perm/regression_SL_lag%d-%d_permPhase%04d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'b','r2','lags','rnames','keptT');
    clear b r2 keptT rp
end
