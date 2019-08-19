clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
iters=1000;

lags=-10:10;
for ei=3;%1:4;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker.mat'],'data','keptvox');
    roi_voxn=1;
    roi_tn=size(gdata,2);
    keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
    
    b=[];
    r2=[];
    F=[];
    for iter=1:iters;
        [data_perm, rp(:,iter)]=phase_rand2(data',1);
        data_perm=data_perm';
        
        for vi=1:size(data,1);
            clear data_mat
            
            %% isc: average across voxels;
            % then zscore overtime, so the results won't be dominated by
            % subjects with higher absolute values.
            
            g=nanmean(gdata(vi,:,:),3);
            y=g;
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=data_perm(vi,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
            
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b(vi,:,iter),~,r(vi,:,iter),~,stats]=regress(y,X);
            
            r2(vi,:,iter)=stats(1);
            F(vi,:,iter)=stats(2);
            p(vi,:,iter)=stats(3);
            
            Y(vi,:,iter)=X*b(vi,:,iter)';
            %  Y=reshape(Y,roi_voxn,length(keptT));
            %   coupling(vi,:)=zeros(1,roi_tn);
            %   coupling(vi,keptT)=corr_col(g(:,keptT),Y);
            
            clear X
        end
    end
    b_sl_perm=b;
    F_sl_perm=F;
    r2_sl_perm=r2;
    p_sl_perm=p;
    r_sl_perm=r;
    lags_sl_perm=lags;
    
    keptT_sl_perm=keptT;
    Y_sl_perm=Y;
    
    % couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm'],'b_sl_perm','F_sl_perm','r2_sl_perm','r_sl_perm','p_sl_perm','lags_sl_perm','keptT_sl_perm','Y_sl_perm','keptvox');
    
    clear b F p r2 rnames coupling r Y
end

toc
