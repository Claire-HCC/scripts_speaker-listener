function roi_tr_bined_pattern_regression_SL_permSubj(perm)

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize_tested=[10 15 30]; % tr;
lags_tested={ -10:-4, -20:-4, -30:-4, -10:-1};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because there is clearly a speech-start effect in the
% listeners' data

for ei= 1:4;
    exp=experiments{ei};
    
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rnames{1} '.mat'],'data');
    tn=size(data,2);
    
    for binSizei=1:length(binSize_tested);
        binSize=binSize_tested(binSizei);
        
        for lagi=1:length(lags_tested);
            lags=lags_tested{lagi};
            
            b=[];
            roi_ids=[];
            r2=nan([length(rnames) tn ]);
            F=nan([length(rnames) tn ]);
            
            for ri=1:length(rnames);
                clear data_mat
                rname=rnames{ri};
                
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
                
                roi_voxn=size(gdata,1);
                tn=size(gdata,2);
                
                data_perm=gdata(:,:,perm);
                gdata_perm=gdata;
                gdata_perm(:,:,perm)=data;
                
                g=nanmean(gdata_perm(:,:,:),3);
                
                for t=1:tn;
                    t_bin=t:(t+binSize-1);
                    
                    if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                        % substract the global mean pattern
                        y=g(:,t_bin);
                        y=y(:);
                        
                        for li=1:length(lags);
                            X(:,:,li)=data_perm(:,t_bin+lags(li));
                        end
                        
                        X=reshape(X,roi_voxn*length(t_bin),length(lags));
                        
                        % centralized X
                        X=X-mean(X);
                        
                        % add an constant
                        
                        [b(ri,t,:),~,~,~,stats]=regress(y,[ones(size(X,1),1) X]);
                        
                        r2(ri,t)=stats(1);
                        F(ri,t)=stats(2);
                        p(ri,t)=stats(3);
                        
                        clear X
                    end
                end
            end
            
            save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/SLg/perm/regression_SL_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' num2str(perm)],'b','F','r2','p','lags','rnames','binSize');
            clear b F p r2 coupling r
        end
    end
end
