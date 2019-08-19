
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
iters=1000;

lags=-4:4;
for ei=1:4;
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
    
    
    for ri=1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/speaker_' rname '.mat'],'data');
        gdata=zscore(nanmean(gdata,1),0,2);
        data=zscore(nanmean(data,1),0,2);
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        g=nanmean(gdata,3);
        y=g;
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
        y=y(:,keptT);
        y=y(:);
        
        for iter=1:iters;
            
            
            %% isc: average across voxels;
            % then zscore overtime, so the results won't be dominated by
            % subjects with higher absolute values.
            
            if ri==1;
                [data_perm, rp(:,iter)]=phase_rand2(data',1);
            else;
                [data_perm]=phase_rand2(data',1,rp(:,iter));
            end
            data_perm=data_perm';
            
            
            for li=1:length(lags);
                X(:,:,li)=data_perm(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b_sl_perm(ri,:,iter),~,r_sl_perm(ri,:,iter),~,stats]=regress(y,X);
            
            r2_sl_perm(ri,:,iter)=stats(1);
            F_sl_perm(ri,:,iter)=stats(2);
            p_sl_perm(ri,:,iter)=stats(3);
            
            Y_sl_perm(ri,:,iter)=X*b_sl_perm(ri,:,iter)';
            %  Y=reshape(Y,roi_voxn,length(keptT));
            %   coupling(ri,:)=zeros(1,roi_tn);
            %   coupling(ri,keptT)=corr_col(g(:,keptT),Y);
            
            
            clear X
        end
    end
    
    
    lags_sl_perm=lags;
    rnames_sl_perm=rnames;
    keptT_sl_perm=keptT;
    
    
    % couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],...
        'b_sl_perm','F_sl_perm','r2_sl_perm','r_sl_perm','p_sl_perm','lags_sl_perm','rnames_sl_perm','keptT_sl_perm','Y_sl_perm','iters');
    clear b_sl_perm F_sl_perm p_sl_perm r2_sl_perm rnames_sl_perm r_sl_perm Y_sl_perm rp
end

toc
