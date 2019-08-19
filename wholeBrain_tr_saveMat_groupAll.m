clear all;

loc='cluster';

set_parameters;
timeUnit='tr' ;

for ei = 2;
    exp=experiments{ei};
    tic % 15 min
     
    if ei<=2; subjn=48; elseif ei>2; subjn=18; end
    
    for s = 1:subjn;
        f=sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listener%02d.nii',expdir,exp,timeUnit,s);
        
        if s==1;
            gdata=nii2mat(f);
        elseif s==10 & strmatch(exp,'pieman'); % the scanning session stopped in the early for this subject
            gdata(:,:,s)=NaN;
        else
            gdata(:,:,s)=nii2mat(f);
        end
    end
    [mat,~]=nii2mat(sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/speaker01.nii',expdir,exp,timeUnit));
    data(:,:,1)=mat;
    
    keptvox=find(sum(sum(gdata==0,3),2)==0 & sum(data==0,2)==0);
    
    roimask=zeros(voxn,1);
    roimask(keptvox)=1;
    save([expdir '/roi_mask/group_fmri_mask_' exp '.mat'],'roimask');
    nii=mat2nii(roimask);
    save_nii(nii,[expdir '/roi_mask/group_fmri_mask_' exp '.nii']);
    
    gdata=gdata(keptvox,:,:);
    save([expdir '/'  exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll.mat'],'gdata','keptvox','-v7.3');
    
    data=data(keptvox,:,:);
    save([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker.mat'],'data','keptvox','-v7.3');
    
    gdata=zscore(gdata,0,2);
    save([expdir '/'  exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_listenerAll.mat'],'gdata','keptvox','-v7.3');
    
    data=zscore(data,0,2);
    save([expdir '/'  exp '/fmri/timeseries/' timeUnit '/wholeBrain/zscore_speaker.mat'],'data','keptvox','-v7.3');
    
   clear data gdata
   toc
end

beep