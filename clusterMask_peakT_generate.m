clear all
loc='mypc'
set_parameters;
exp='merlin';
subj='listenerZscoreMean_g2';
relation='LL';
p_thr=0.01;
rname='HG_L';


lags=-10:10;

peakT_nii=load_nii([expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '/isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags)) '_' rname '_peakT_p' num2str(p_thr) 'FDR.nii']);

labels=zeros(volsize);
for lagi=1:length(lags);
    lag=lags(lagi);
    img_temp=(peakT_nii.img);
    
    % label cluster
    label=zeros(volsize);
    img_temp_bw=(img_temp==lag);
    
    % cluster_properties = regionprops(img_temp_bw,'area','pixellist');
    cluster_properties = bwconncomp(img_temp_bw);
    % cluster_size=[cluster_properties(:).Area];
    cluster_size=regionprops(cluster_properties ,'Area');
    cluster_size=[cluster_size(:).Area];
    smallClusterI=find(cluster_size<30);
    for ci=1:length(smallClusterI);
        sci=smallClusterI(ci);
        img_temp_bw(cluster_properties.PixelIdxList{sci})=0;
    end
    
    label=bwlabeln(img_temp_bw,26);
    % plus 0.001 for lag 0
    label(label>0)=(abs(lag)*100+label(label>0))*sign(lag+0.001);
    labels=nansum(cat(4,label,labels),4);
end

labels_nii=peakT_nii;
labels_nii.img=labels;
save_nii(labels_nii,[expdir exp '/fmri/nii/wholeBrain/isfc_seed/' relation '//isfc_seed_lagcorr_'  num2str(min(lags)) '-' num2str(max(lags)) '_' rname '_peakT_p' num2str(p_thr) 'FDR_clusterMask.nii'])