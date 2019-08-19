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
    keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));
    
    b=[];
    r2=[];
    F=[];
    for vi=1:size(data,1);
        clear data_mat
        
        for s=1:size(gdata,3);
            
            othersi=1:size(gdata,3);
            others=othersi(othersi~=s);
            others=nanmean(gdata(vi,:,othersi),3);
            
            self=gdata(vi,:,s);
            
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
            
            [b(vi,:,s),~,r(vi,:,s),~,stats]=regress(y,X);
            
            r2(vi,:,s)=stats(1);
            F(vi,:,s)=stats(2);
            p(vi,:,s)=stats(3);
            Y(vi,:,s)=X*b(vi,:,s)';
            clear X
        end
        
    end
    
    b_ll=b;
    F_ll=F;
    r2_ll=r2;
    p_ll=p;
    r_ll=r;
    lags_ll=lags;
    keptT_ll=keptT;
    Y_ll=Y;
    
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/' '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags))],'b_ll','F_ll','r2_ll','p_ll','r_ll','lags_ll','keptT_ll','Y_ll','keptvox','-v7.3');
    
    mat=zeros(voxn,1);
    mat(keptvox)=nanmean(F_ll,3);
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_F.nii']);
    
    mat=zeros(voxn,1);
    mat(keptvox)=nanmean(r2_ll,3);
    nii=mat2nii(mat);
    save_nii(nii,[expdir '/' exp '/fmri/isc_regression/' timeUnit '/wholeBrain/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_r2.nii']);

    
    clear b F p r2  coupling r Y
end

toc
