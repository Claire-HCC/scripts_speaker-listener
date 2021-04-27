function fc_cluster
% use cluster from pieman old to replace mor.
% i probably will have to define my own roi...
loc='mypc';
set_parameters;
timeUnit='tr' ;
k=6;
perc=0.30;
kmeansd='cityblock';
iscOverlapPercThr=0.75; % the percentage of exepriments showing strong isc

load([expdir '/pieman_rest/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc.mat'],'rzm','keptvox','keptT');

% masked by story isc
f= sprintf('%s/roi_mask/isc%dPercMasks.mat',expdir,perc*100);
load(f,'roimask');
roimask=roimask>=iscOverlapPercThr;
rzm=rzm(ismember(keptvox,find(roimask)),ismember(keptvox,find(roimask)));
keptvox=keptvox(ismember(keptvox,find(roimask)));

[idx, C,sumd, D]= kmeans(tanh(rzm),k,'Distance',kmeansd);

save([expdir '/pieman_rest/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(iscOverlapPercThr*100) 'Overlap_cluster' num2str(k) '.mat'],'C','idx','D','keptvox','keptT','k','kmeansd','-v7.3');

mat=zeros(voxn,1);
mat(keptvox)=idx;
nii=mat2nii(mat);
save_nii(nii,[expdir '/pieman_rest/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k) '.nii'])

figure;
[~,sorti]=sort(idx);
imagesc(rzm(sorti,sorti),[-0.2 0.3]);


%% check network manually here. Then set the network order
% networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
load([expdir '/' 'pieman_rest' '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k) '.mat' ],'C','idx','D','keptvox','keptT','k');
idx_newOrder=[1 2 5 6 4 3];
C=C(idx_newOrder,:);
D=D(:,idx_newOrder);
idx_new=zeros(size(idx));
for ki=1:k;
    idx_new(idx==idx_newOrder(ki))=ki;
end
idx=idx_new;
mat=zeros(voxn,1);
mat(keptvox)=idx;
nii=mat2nii(mat);
save_nii(nii,[expdir '/' 'pieman_rest' '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k)  '.nii']);
save([expdir '/' 'pieman_rest' '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k) ],'C','idx','D','keptvox','keptT','k','networks');


%% create masks
mkdir([expdir '/roi_mask/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k) '/nii/']);
mkdir([expdir '/roi_mask/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k) '/mat/']);
for ki=1:k;
    roimask=zeros(voxn,1);
    roimask(keptvox(idx==ki))=1;
    
    save([expdir '/roi_mask/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k) '/mat/' networks{ki} '.mat'],'roimask');
    
    nii=mat2nii(roimask);
    save_nii(nii,[expdir '/roi_mask/restFc_isc' num2str(100*perc) 'PercMasked_' num2str(100*iscOverlapPercThr) 'Overlap_cluster' num2str(k)  '/nii/' networks{ki} '.nii']);
end
