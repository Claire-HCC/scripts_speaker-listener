clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
binSize=30; % tr;
binStep=1;
lags=-60:-4;

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

for ei=3;%1:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    roi_ids=[];
    
    b_l=[];
    r2_l=[];
    F_l=[];
    b_ls=[];
    r2_ls=[];
    F_ls=[];
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0])-binSize+1);
        binN=length(1:binStep:(length(keptT)-binSize+1));
        
        g=nanmean(gdata,3);
        
        for bi=1:binN;
            keptTi=((bi-1)*binStep+1):((bi-1)*binStep+binSize);
            
            y=g;
            y=y(:,keptT(keptTi));
            y=y(:);
            
            for li=1:length(lags);
                Xl(:,:,li)=g(:,keptT(keptTi)+lags(li));
                Xs(:,:,li)=data(:,keptT(keptTi)+lags(li));
            end
            
            Xl=reshape(Xl,roi_voxn*length(keptTi),length(lags));
            Xs=reshape(Xs,roi_voxn*length(keptTi),length(lags));
            
            % centralized X
            Xl=Xl-mean(Xl);
            Xs=Xs-mean(Xs);
            
            % add an constant
            Xl=[ones(size(Xl,1),1) Xl];
            
            [b_l(ri,keptT(keptTi(1)),:),~,~,~,stats_l]=regress(y,Xl);
            [b_ls(ri,keptT(keptTi(1)),:),~,~,~,stats_ls]=regress(y,[Xl Xs]);
            
            r2_l(ri,keptT(keptTi(1)))=stats_l(1);
            F_l(ri,keptT(keptTi(1)))=stats_l(2);
            p_l(ri,keptT(keptTi(1)))=stats_l(3);
            
            r2_ls(ri,keptT(keptTi(1)))=stats_ls(1);
            F_ls(ri,keptT(keptTi(1)))=stats_ls(2);
            p_ls(ri,keptT(keptTi(1)))=stats_ls(3);
            
            
            clear Xl Xs
        end
    end
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_l','b_ls','F_l','r2_l','p_ls','F_ls','r2_ls','p_ls','keptT','rnames');
    % clear b F p r2 rnames coupling
end

toc
beep
