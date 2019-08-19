clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;

lags=-10:10;
for ei=1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat' ],'gdata');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker.mat'],'data','keptvox');
    roi_voxn=1;
    roi_tn=size(gdata,2);
            keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
            
    b=[];
    r2=[];
    F=[];
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
            X(:,:,li)=data(vi,keptT+lags(li));
        end
        
        X=reshape(X,roi_voxn*length(keptT),length(lags));
        
        % centralized X
        X=X-mean(X);
        
        % add an constant
        X=[ones(size(X,1),1) X];
        
        [b(vi,:),~,r(vi,:),~,stats]=regress(y,X);
        
        r2(vi,:)=stats(1);
        F(vi,:)=stats(2);
        p(vi,:)=stats(3);
        
        Y(vi,:)=X*b(vi,:)';
        %  Y=reshape(Y,roi_voxn,length(keptT));
        %   coupling(vi,:)=zeros(1,roi_tn);
        %   coupling(vi,keptT)=corr_col(g(:,keptT),Y);
        
        clear X
    end
    
    b_sl=b;
    F_sl=F;
    r2_sl=r2;
    p_sl=p;
    r_sl=r;
    lags_sl=lags;
    
    keptT_sl=keptT;
    Y_sl=Y;
    
    % couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_sl','F_sl','r2_sl','r_sl','p_sl','lags_sl','keptT_sl','Y_sl','keptvox');
    mat=zeros(voxn,1);
    mat(keptvox)=F_sl;
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_F.nii']);
    
       mat=zeros(voxn,1);
    mat(keptvox)=r2_sl;
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_r2.nii']);

    
    clear b F p r2 rnames coupling r Y
end

toc
