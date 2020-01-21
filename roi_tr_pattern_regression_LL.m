clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
lags_tested={-10:10, -20:20, -30:30, 0:0};

for ei=2:4;
    exp=experiments{ei};
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for lagi=1:length(lags_tested);
        lags=lags_tested{lagi};
        
        for ri=1:length(rnames);
            clear data_mat
            rname=rnames{ri};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
            
            roi_voxn=size(gdata,1);
            roi_tn=size(gdata,2);
            keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));
            
            for s=1:size(gdata,3)
                othersi=1:size(gdata,3);
                othersi=othersi(othersi~=s);
                others=nanmean(gdata(:,:,othersi),3);
                
                self=gdata(:,:,s);
                
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
                r2_byTime(ri,:,s)=nan(roi_tn,1);
                r2_byTime(ri,keptT,s)=1-ssr./sst;
                
                clear X
                
            end
        end
        
        save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT','r2_byTime');
        clear b F p r2 r2_byTime
    end
    
end

toc
