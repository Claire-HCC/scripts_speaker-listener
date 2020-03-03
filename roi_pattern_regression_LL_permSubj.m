function roi_tr_pattern_regression_LL_permSubj(perm);

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
lags=0;

for ei=2;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rnames{1} '.mat' ],'gdata')
    [~ ,tn ,subjn]=size(gdata);
    
    if ismember(perm,find(~isnan(gdata(1,1,:))));
        b=nan([length(rnames) length(lags)+1 subjn]);
        r2=nan([length(rnames)  subjn]);
        F=nan([length(rnames)  subjn]);
        p=nan([length(rnames) subjn]);
        r2_byTime=nan([length(rnames),tn, subjn]);
        
        for ri=1:length(rnames);
            clear data_mat
            rname=rnames{ri};
            
            if  exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat' ],'data');
                
                roi_voxn=size(gdata,1);
                keptT=(abs(min(lags))+1+crop_start):(tn-max([lags,0]));
                
                gdata_perm=gdata;
                gdata_perm(:,:,perm)=data;
                
                for s=1:size(gdata_perm,3)
                    othersi=1:size(gdata_perm,3);
                    othersi=othersi(othersi~=s);
                    others=nanmean(gdata_perm(:,:,othersi),3);
                    
                    self=gdata_perm(:,:,s);
                    
                    y=self;
                    y=y(:,keptT);
                    y=y(:);
                    
                    for li=1:length(lags);
                        X(:,:,li)=others(:,keptT+lags(li));
                    end
                    
                    X=reshape(X,roi_voxn*length(keptT),length(lags));
                    
                    % centralized X
                    X=X-mean(X);
                    
                    [b(ri,:,s),~,r,~,stats]=regress(y,[ones(size(X,1),1) X]);
                    
                    r2(ri,s)=stats(1);
                    F(ri,s)=stats(2);
                    p(ri,s)=stats(3);
                    
                    ymat=reshape(y,roi_voxn,length(keptT));
                    rmat=reshape(r,roi_voxn,length(keptT));
                    
                    ssr=sum(rmat.^2);
                    sst=sum((ymat-mean(y)).^2);
                    r2_byTime(ri,:,s)=nan(tn,1);
                    r2_byTime(ri,keptT,s)=1-ssr./sst;
                    
                    clear X
                end
            end
        end
        
        save(sprintf('%s/%s/fmri/pattern_regression/%s/roi/%s/LLselfother/perm/regression_LL_lag%d-%d_permSL%03d',expdir,exp,timeUnit,froidir,min(lags),max(lags),perm),'b','F','r2','p','lags','rnames','keptT','r2_byTime');
        clear b F p r2 rnames r2_byTime
    end
end

