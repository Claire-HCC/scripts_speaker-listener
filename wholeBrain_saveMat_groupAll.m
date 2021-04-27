clear all;

loc='cluster';
set_parameters;
timeUnit='tr_uncropped' ;
load([expdir 'roi_mask/MNI152NLin2009cAsym_3x3x4mm_brain_gm.mat'],'roimask')
mask_gray=roimask;

for ei =1;
    exp=exp_parameters.experiments{ei};
    tic % 15 min
    
    fs=dir(sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/*.nii',expdir,exp,timeUnit));
    fs={fs.name};
    
    subjs=sscanf(sprintf('%s',fs{:}),'listener%d.nii')
    [~,subjsOrder]=sort(subjs)
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/%s.nii',expdir,exp,timeUnit,fs{1});
    gdata=nii2mat(f);
    
    keptvox=find(mask_gray==1);
    
    gdata=nan(length(keptvox),size(gdata,2),length(subjs));
    
    for si=1:length(subjs);
        %    s=sscanf(fs{si}, 'sub-%d');
        s=subjsOrder(si);
        f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/%s.nii',expdir,exp,timeUnit,fs{s});
        mat=nii2mat(f);
        gdata(:,:,si)=mat(keptvox,:);
    end
    
    gmask=find(~isnan(nanmean(nanmean(gdata,3),2)) & nanmean(nanmean(gdata,3),2)~=0);
    gdata=gdata(gmask,:,:);
    keptvox=keptvox(gmask);

    save([expdir '/'  exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll.mat'],'gdata','fs','keptvox','-v7.3');
    
%     gdata=zscore(gdata,0,2);
%     save([expdir '/'  exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_zscore.mat'],'gdata','fs','keptvox','-v7.3');
%     
%     roimask=zeros(voxn,1);
%     roimask(keptvox)=1;
%     save([expdir '/roi_mask/group_fmri_mask_' exp '.mat'],'roimask');
%     nii=mat2nii(roimask);
%     save_nii(nii,[expdir '/roi_mask/group_fmri_mask_' exp '.nii']);
%     
    clear gdata
end
toc


%
% % speaker
% for ei = 3:8;% 1:2;%1:4;
%     exp=experiments{ei};
%     tic % 15 min
%
%     [mat,~]=nii2mat(sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/speaker01.nii',expdir,exp,timeUnit));
%     data(:,:,1)=mat;
%
%     keptvox=find(sum(sum(gdata==0,3),2)==0 & mask_gray==1);
%
%     roimask=zeros(voxn,1);
%     roimask(keptvox)=1;
%     save([expdir '/roi_mask/group_fmri_mask_' exp '.mat'],'roimask');
%     nii=mat2nii(roimask);
%     save_nii(nii,[expdir '/roi_mask/group_fmri_mask_' exp '.nii']);
%
%     data=data(keptvox,:,:);
%     save([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker.mat'],'data','keptvox','-v7.3');
%
%     data=zscore(data,0,2);
%     save([expdir '/'  exp '/fmri/timeseries/' timeUnit '/wholeBrain/speaker_zscore.mat'],'data','keptvox','-v7.3');
%
%     clear data gdata
%     toc
% end
%
% beep
% beep