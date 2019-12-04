clear all;
% 0.5 hr for 1 roi 1000 permutation
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
lags=0;
iters=1000;

for ei=3;%[1 2 4]%1:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    
    
    ris=find(ismember(rnames,{'vPCUN'}))
    for i=1:length(ris);%1:length(rnames);
        ri=ris(i);
        rname=rnames{ri};
        
        clear data_mat
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        eventLabels_LG=h5read([expdir   exp '/fmri/hmm/' rname  '_findListenersEventInSpeaker.hdf5'],'///eventLabels_LG');
        K=length(unique(eventLabels_LG(:,1)));
        eventLabels_LG(1:(1-min(lags)))=0;
        eventLabels_LG_kept=unique(eventLabels_LG(eventLabels_LG~=0));
        
        
        for iter=1:iters;
            for s=1:size(gdata,3);
                if i==1 & s==1;
                    [gdata_perm, rp(:,iter)]=phase_rand2(gdata(:,:,s)',1);
                    gdata_perm=gdata_perm';
                else;
                    gdata_perm(:,:,s)=phase_rand2(gdata(:,:,s)',1,rp(:,iter))';
                end
            end
            
            
            listeners=1:size(gdata_perm,3);
            for s=listeners;
                
                othersi=listeners(listeners~=s);
                others=nanmean(gdata_perm(:,:,othersi),3);
                self=gdata_perm(:,:,s);
                
                for evi=1:K;
                    keptT=find(eventLabels_LG==evi);
                    
                    y=others;
                    y=y(:,keptT);
                    y=y(:);
                    
                    for li=1:length(lags);
                        X(:,:,li)=self(:,keptT+lags(li));
                    end
                    
                    X=reshape(X,roi_voxn*length(keptT),length(lags));
                    
                    % centralized X
                    X=X-mean(X);
                    
                    % add an constant
                    X=[ones(size(X,1),1) X];
                    
                    [b_perm(ri,:,s,iter),~,~,~,stats]=regress(y,X);
                    
                    r2_perm(ri,evi,s,iter)=stats(1);
                    F_perm(ri,evi,s,iter)=stats(2);
                    p_perm(ri,evi,s,iter)=stats(3);
                    
                    clear X
                end
                
            end
        end
    end
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'b_perm','F_perm','r2_perm','p_perm','lags','eventLabels_LG_kept','rnames','-v7.3');
    clear b_perm F_perm p_perm r2_perm rnames
end

toc
