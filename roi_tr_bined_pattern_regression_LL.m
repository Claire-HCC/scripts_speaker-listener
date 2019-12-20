function roi_tr_bined_pattern_regression_LL

loc='cluster';
tic
set_parameters;
timeUnit='tr' ;

froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

% cropt start because ther eis clearly a spech-start effect in the
% listeners' data
binSize_tested=[10 15 30]; % tr;
lags=0;

for binSizei=1:length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    for ei=[1 2 4];
        exp=experiments{ei};
        rnames=table2array(roi_table(:,3));
        ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
        rnames=rnames(ris);
        
        for ri=1:length(rnames);
            clear data_mat
            rname=rnames{ri};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
            roi_voxn=size(gdata,1);
            tn=size(gdata,2);
            listeners=1:size(gdata,3);
            
            for s=listeners;
                othersi=listeners(listeners~=s);
                others=nanmean(gdata(:,:,othersi),3);
                self=gdata(:,:,s);
                
                for t=1:tn;
                    t_bin=t:(t+binSize-1);
                    
                    if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                        
                        y=others;
                        y=y(:,t_bin);
                        y=y(:);
                        
                        for li=1:length(lags);
                            X(:,:,li)=self(:,t_bin+lags(li));
                        end
                        
                        X=reshape(X,roi_voxn*length(t_bin),length(lags));
                        
                        % centralized X
                        X=X-mean(X);
                        
                        [b(ri,t,:,s),~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
                        
                        r2(ri,t,s)=stats(1);
                        F(ri,t,s)=stats(2);
                        p(ri,t,s)=stats(3);
                        
                    else
                        b(ri,t,s)=nan;
                        r2(ri,t,s)=nan;
                        F(ri,t,s)=nan;
                        p(ri,t,s)=nan;
                    end
                    clear X y
                end
            end
            disp(ri)
        end
        save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','-v7.3');
        clear b F r2 p rnames
    end
end
toc





