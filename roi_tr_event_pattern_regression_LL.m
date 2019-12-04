clear all;

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
        
        listeners=1:size(gdata,3);
        for s=listeners;
            
            othersi=listeners(listeners~=s);
            others=nanmean(gdata(:,:,othersi),3);
            self=gdata(:,:,s);
            
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
                
                [b(ri,:,s),~,~,~,stats]=regress(y,X);
                
                r2(ri,evi,s)=stats(1);
                F(ri,evi,s)=stats(2);
                p(ri,evi,s)=stats(3);
                
                Y=X*b(ri,:,s)';
                Y=reshape(Y,roi_voxn,length(keptT));
                y=reshape(y,roi_voxn,length(keptT));
                
                clear X
            end
            
        end
    end
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','eventLabels_LG_kept','y','Y','rnames','-v7.3');
    clear b F p r2 rnames coupling
end

toc
