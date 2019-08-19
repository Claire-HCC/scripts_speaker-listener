
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

lags=-5:5;
for ei=1%:2;
    exp=experiments{ei};
    b_SL=[];
    F=[];
    r2=[];
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain//zscore_listenerAll'  ],'gdata','keptvox');
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain//zscore_speaker'  ],'data','keptvox');
    
    roi_voxn=size(gdata,1);
    roi_tn=size(gdata,2);
    
    g=nanmean(gdata,3);
    keptT=(max(lags)+1):(size(data,2)-max(lags));
    for vi=1:length(keptvox);
        
        y=g(vi,keptT);
        y=y(:);
        
        for li=1:length(lags);
            X(:,:,li)=data(vi,keptT+lags(li));
        end
        
        X=reshape(X,length(keptT),length(lags));
        X=[ones(size(X,1),1) X];
        
        [b(vi,:),~,~,~,stats]=regress(y,X);
        
        r2(vi,:)=stats(1);
        F(vi,:)=stats(2);
        p(vi,:)=stats(3);
        
        Y=X*b(vi,:)';
        Y=reshape(Y,1,length(keptT));
        
        coupling(vi,:)=zeros(1,roi_tn);
        coupling(vi,keptT)=corr_col(g(vi,keptT),Y);
        
        clear X
        
    end
end
 couplingz=0.5*log((1+coupling)./(1-coupling));
save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/wholeBrain/regression_SL' ],'b','F','r2','p','couplingz','lags','keptT','keptvox');

mat=zeros(voxn,1);
mat(keptvox)=F;
nii=mat2nii(mat);
save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/wholeBrain/regression_SL_F.nii' ])
% clear b_SL roi_ids rnames
mat=zeros(voxn,1);
mat(keptvox)=r2;
nii=mat2nii(mat);
save_nii(nii,[expdir '/' exp '/fmri/pattern_regression/' timeUnit '/wholeBrain/regression_SL_r2.nii' ])

toc
beep