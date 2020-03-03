function roi_tr_temporal_regression_SL(ri)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-4:4};

for ei=1:4;%2:4;
    exp=experiments{ei};
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rnames{1} '.mat' ],'data')
        tn=size(data,2);
        
        b=nan([length(rnames) length(lags)+1 ]);
        r2=nan([length(rnames) 1 ]);
        F=nan([length(rnames) 1]);
        p=nan([length(rnames) 1]);
        r2_byTime=nan([length(rnames),tn]);
        
            rname=rnames{ri};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/permPhase/speaker_' rname '.mat'],'data');
                data_perm=data;
                roi_voxn=size(gdata,1);
                tn=size(data,2);
                
                keptT_s=find(([1:tn]+min(lags))==1)+crop_start;
                keptT_e=min(tn,tn-max(lags));
                keptT=keptT_s:keptT_e;
                
                gdata=nanmean(gdata,1);
                g=nanmean(gdata,3);
                y=g;
                y=y(:,keptT);
                y=y(:);
                
                for perm=1:permN;
                    data=data_perm(:,:,perm);
                    
                    % average all voxels within the roi
                    data=nanmean(data,1);
                    
                    for li=1:length(lags);
                        X(:,:,li)=data(:,keptT+lags(li));
                    end
                    
                    X=reshape(X,roi_voxn*length(keptT),length(lags));
                    X=X-mean(X);
                    
                    [b(ri,:),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                    
                    r2(ri,1)=stats(1);
                    F(ri,1)=stats(2);
                    p(ri,1)=stats(3);
                    
                    ymat=reshape(y,roi_voxn,length(keptT));
                    rmat=reshape(r,roi_voxn,length(keptT));
                    ssr=sum(rmat.^2);
                    sst=sum((ymat-mean(y)).^2);
                    
                    r2_byTime(ri,:)=nan(tn,1);
                    r2_byTime(ri,keptT)=1-ssr./sst;
                    
                    clear X
                end
            end
        
        save([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLg/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_' rname],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
        clear b F p r2 r2_byTime
    end
end

