
clear all;

loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

tic % 15 min
iters=1000;

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;

lags=-10:10;
for ei=1:2;%1:4;
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
        
        rname=rnames{ri};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat'],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
        
        for iter=1:iters;
            clear data_mat
            data_perm=zeros(size(data))';
            if ri==1;
                [data_perm((crop_start+1):end,:), rp(:,iter)]=phase_rand2(data(:,(crop_start+1):end)',1);
            else;
                [data_perm((crop_start+1):end,:)]=phase_rand2(data(:,(crop_start+1):end)',1,rp(:,iter));
            end
            data_perm=data_perm';
            
            g=nanmean(gdata,3);
            y=g;
            
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data_perm(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b_perm(ri,:,iter),~,~,~,stats]=regress(y,X);
            r2_perm(ri,:,iter)=stats(1);
            F_perm(ri,:,iter)=stats(2);
            p_perm(ri,:,iter)=stats(3);
            
            Y=X*b_perm(ri,:,iter)';
            Y=reshape(Y,roi_voxn,length(keptT));
            coupling_perm(ri,:,iter)=zeros(1,roi_tn);
            coupling_perm(ri,keptT,iter)=corr_col(g(:,keptT),Y);
            
            clear X
            
        end
        
    end
    couplingz_perm=0.5*log((1+coupling_perm)./(1-coupling_perm));
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'couplingz_perm','b_perm','F_perm','r2_perm','p_perm','lags','rnames');
    clear p_perm b_perm F_perm r2_perm rnames coupling_perm rp
end

toc
