
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
crop_start=10;

lags=-10:10;
for ei=1%:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    for ri=3%1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat'],'data');
        
        %% isc: average across voxels; 
        % then zscore overtime, so the results won't be dominated by
        % subjects with higher absolute values.
        gdata=zscore(nanmean(gdata,1),0,2);
        data=zscore(nanmean(data,1),0,2);
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        g=nanmean(gdata,3);
        y=g;
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
        y=y(:,keptT);
        y=y(:);
        
        for li=1:length(lags);
            X(:,:,li)=data(:,keptT+lags(li));
        end
        
        X=reshape(X,roi_voxn*length(keptT),length(lags));
        
        % centralized X
        X=X-mean(X);
        
        % add an constant
        X=[ones(size(X,1),1) X];
        
        [b(ri,:),~,r(ri,:),~,stats]=regress(y,X);
        
        r2(ri,:)=stats(1);
        F(ri,:)=stats(2);
        p(ri,:)=stats(3);
        
        Y(ri,:)=X*b(ri,:)';
        %  Y=reshape(Y,roi_voxn,length(keptT));
        %   coupling(ri,:)=zeros(1,roi_tn);
        %   coupling(ri,keptT)=corr_col(g(:,keptT),Y);
        
        
        clear X
    end
    
    b_sl=b;
    F_sl=F;
    r2_sl=r2;
    p_sl=p;
    r_sl=r;
    lags_sl=lags;
    rnames_sl=rnames;
    keptT_sl=keptT;
    Y_sl=Y;
    
    % couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_sl','F_sl','r2_sl','r_sl','p_sl','lags_sl','rnames_sl','keptT_sl','Y_sl');
    clear b F p r2 rnames coupling r Y
end

toc
