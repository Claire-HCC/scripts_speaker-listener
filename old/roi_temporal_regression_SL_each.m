function roi_temporal_regression_SL_each(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));

crop_start=10;
lags_tested={-40:40,-4:4,-10:-4. 4:10, -10:10, -10:-1, 0, 1:10};

exp=experiments{ei};
% mkdir([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/perm']);

for lagi=1;%1:length(lags_tested);
    lags=lags_tested{lagi};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rnames{1} '.mat' ],'gdata')
    [~,tn,listenerN]=size(gdata);
    
    b=nan([length(rnames) length(lags) listenerN]);
    r2=nan([length(rnames) 1 listenerN]);
    F=nan([length(rnames) 1 listenerN]);
    p=nan([length(rnames) 1 listenerN]);
    r2_byTime=nan([length(rnames),tn, listenerN]);
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_zscore_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_zscore_' rname '.mat'],'data');
            % gdata(:,:,subjects_excluded{ei})=NaN;
            data=nanmean(data,1);
            gdata=nanmean(gdata,1);
            
            roi_voxn=1;
            tn=size(data,2);
            
            keptT_s= max(-min(lags),0)+crop_start+1;
            keptT_e=min(tn,tn-max(lags));
            keptT=keptT_s:keptT_e;
            
            for li=1:length(lags);
                X(:,:,li)=data(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            X=X-mean(X);
            
            for si=1:listenerN;
                y=gdata(:,keptT,si);
                y=y(:);
                
                [b_temp,~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                
                b(ri,:,si)=b_temp(2:end);
                r2(ri,1,si)=stats(1);
                F(ri,1,si)=stats(2);
                p(ri,1,si)=stats(3);
                
                ymat=reshape(y,roi_voxn,length(keptT));
                rmat=reshape(r,roi_voxn,length(keptT));
                ssr=sum(rmat.^2);
                sst=sum((ymat-mean(y)).^2);
                
                r2_byTime(ri,:,si)=nan(tn,1);
                r2_byTime(ri,keptT,si)=1-ssr./sst;
            end
            clear X
        end
    end
    
    save([expdir '/' exp '/fmri/temporal/regression/' timeUnit '/roi/' froidir '/SL_each/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
    clear b F p r2 r2_byTime
end


