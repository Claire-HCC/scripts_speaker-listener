
clear all
close all
set_parameters
thr=0.18;

intervals=[-20 -2;
    -1 1;
    2 20
    ]
interval_names={'l_precede','sl_syn','s_precede'}
for ei=1%:2;
    exp= experiments{ei};
    
    load([expdir exp '/fmri/mat/wholeBrain/' exp '_isc_listener_peakR.mat' ],'isc_peakR');
    load([expdir exp '/fmri/mat/wholeBrain/' exp '_isc_listener_peakT.mat' ],'isc_peakT');
    
    r=zeros(size(isc_peakR));
    voxn=[];
    
    for intvi = 1:size(intervals,1);
        for si=1:size(isc_peakR,3);
            i=find(isc_peakR(:,:,si)>thr & isc_peakT(:,:,si) >= intervals(intvi,1) & isc_peakT(:,:,si) <= intervals(intvi,2));
            r(i,:,si)=isc_peakR(i,:,si);
            voxn(si,intvi)=sum(r(:,:,si)~=0);
        end
        
%         nii.hdr=load_nii_hdr([expdir 'roi_mask/MNI152_T1_3mm_brain.nii']);
%         data=mean(r,3);
%         data((end+1):(61*73*61),:)=0;
%         nii.img=reshape(data,61,73,61);
%         save_nii(nii,[expdir exp '/fmri/nii/wholeBrain/' exp '_' interval_names{intvi} '.nii']);
    clear r
    end
end

load([expdir exp '\bhv_merlin.mat'])
[~,~,rk1]  = unique(bhv.scores);
figure;
 for intvi = 1:size(intervals,1);

[~,~, rk2]  = unique(voxn(:,intvi));
[r p]=corr(bhv.scores,voxn(:,intvi),'type','spearman');
subplot(1,3,intvi);

scatter(rk1,rk2);
title({interval_names{intvi},['Spearman r=' num2str(r) '; p=' num2str(p)]});
 end
 

