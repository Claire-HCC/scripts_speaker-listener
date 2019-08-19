% SL + LL, 1000 combinations takes 7 hr on scotty
clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
bootsN=1000;

lags=-4:4;
for ei=3%1:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bootstrap' ],'g1');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rnames{1} '.mat' ],'gdata');
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    for  booti=1:bootsN;
        subjs=1:size(gdata,3);
        g2=subjs(~ismember(1:size(gdata,3),g1(booti,:)));
        
        for ri=1:length(rnames);
            clear data_mat
            rname=rnames{ri};
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
            
            %% isc or pattern
            % gdata=nanmean(gdata,1);
            
            roi_voxn=size(gdata,1);
            roi_tn=size(gdata,2);
            keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));
            
            for s=g2;
                
                othersi=g2(g2~=s);
                others=nanmean(gdata(:,:,othersi),3);
                
                self=gdata(:,:,s);
                
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
                
                [b(ri,:,s,booti),~,~,~,stats]=regress(y,X);
                
                r2(ri,:,s,booti)=stats(1);
                F(ri,:,s,booti)=stats(2);
                p(ri,:,s,booti)=stats(3);
                
                Y=X*b(ri,:,s,booti)';
                Y=reshape(Y,roi_voxn,length(keptT));
                coupling(ri,:,s,booti)=zeros(1,roi_tn);
                coupling(ri,keptT,s,booti)=corr_col(others(:,keptT),Y);
                
                clear X
            end
        end
    end
    couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bootstrap2'],'b','F','r2','p','couplingz','lags','rnames','keptT','-v7.3');
    clear b F p r2 rnames coupling
end

toc

  save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bootstrap_cp'],'couplingz''-v7.3');
  
