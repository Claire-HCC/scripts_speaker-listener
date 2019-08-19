clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because there is clearly a speech-start effect in the
% listeners' data
crop_start=10;
binSize=30; % tr;
binStep=1;

lags=-10:10;
for ei=1:3;%1:4;
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
    for ri=1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        g=nanmean(gdata(:,:,:),3);
        
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
        binN=length(1:binStep:(length(keptT)-binSize+1));
        
        for bi=1:binN;
            keptTi=((bi-1)*binStep+1):((bi-1)*binStep+binSize);
            
            % substract the global mean pattern
            y=g(:,keptT(keptTi));
            y=y(:);
            for li=1:length(lags);
                X(:,:,li)=data(:,keptT(keptTi)+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptTi),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b(ri,keptT(keptTi(1)),:),~,~,~,stats]=regress(y,X);
            
            r2(ri,keptT(keptTi(1)))=stats(1);
            F(ri,keptT(keptTi(1)))=stats(2);
            p(ri,keptT(keptTi(1)))=stats(3);
            
            %             Y=X*squeeze(b(ri,bi,:));
            %             Y=reshape(Y,roi_voxn,length(keptTi));
            %             y=reshape(y,roi_voxn,length(keptTi));
            %             %  coupling(ri,:)=zeros(1,length(keptT));
            %             coupling(ri,keptT(keptTi))=corr_col(y,Y);
            %
                         clear X
        end
    end
   % couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_SL_bin' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b','F','r2','p','lags','rnames','keptT','binSize','binN','binStep');
    clear b F p r2 rnames coupling r
end

toc
beep
