function fc_cluster
% use cluster from pieman old to replace mor.
% i probably will have to define my own roi...
loc='mypc';
set_parameters;
timeUnit='tr' ;
exp='pieman_rest';
ei=find(ismember(exp_parameters.experiments,exp));
k=6;
perc=0.15;
kmeansd='cityblock';

load([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked'],'rzm','keptvox','keptT');

[idx, C,sumd, D]= kmeans(tanh(rzm),k,'Distance',kmeansd);

save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k)],'C','idx','D','keptvox','keptT','k','kmeansd','-v7.3');

mat=zeros(voxn,1);
mat(keptvox)=idx;
nii=mat2nii(mat);
save_nii(nii,[expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k) '.nii'])

figure;
[~,sorti]=sort(idx);
imagesc(rzm(sorti,sorti),[-0.2 0.3]);


%% check network manually here. Then set the network order
exp='pieman_rest';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
load([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k)],'C','idx','D','keptvox','keptT','k');
idx_newOrder=[1 2 4 3 5 6];
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
save_nii(nii,[expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k) '.nii'])

save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k)],'C','idx','D','keptvox','keptT','k','networks');


%% create masks
mkdir([expdir '/roi_mask/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k) '/nii/']);
mkdir([expdir '/roi_mask/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k) '/mat/']);
for ki=1:k;
    roimask=zeros(voxn,1);
    roimask(keptvox(idx==ki))=1;
    
    save([expdir '/roi_mask/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k) '/mat/network' num2str(ki) '.mat'],'roimask');
    
    nii=mat2nii(roimask);
    save_nii(nii,[expdir '/roi_mask/fc_isc' num2str(100*perc) 'PercMasked_cluster' num2str(k) '/nii/network' num2str(ki) '.nii']);
end

%% use trained clusters
% for eii=2;%[7];%1:length(eis);
%     ei=eis(eii);
%     exp=experiments{ei};
%
%     load([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc'],'rzm','keptvox','keptT');
%     keptvox_test=keptvox;
%     keptvox=intersect(keptvox_ref,keptvox_test);
%
%     vi=find(ismember(keptvox_test,keptvox));
%     [~,idx] = pdist2(C_ref(:,ismember(keptvox_ref,keptvox)),rzm(vi,vi),kmeansd,'Smallest',1);
%     [~,sorti]=sort(idx);
%
%     imagesc(rzm(vi(sorti),vi(sorti)))
%     mat=zeros(voxn,1);
%     mat(keptvox)=idx;
%     nii=mat2nii(mat);
%     save_nii(nii,'test.nii')
%     %    save_nii(nii,[expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_cluster.nii'])
%
%     %  save([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_cluster'],'C','idx','D','keptvox','keptT','-v7.3');
% end

%% compute dice index
% % the first one as the reference
% eis=[1 2  11 12 9 10 ];
% load([expdir '/pieman_old/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_cluster'],'C','idx','keptvox','keptT');
% idx_ref=idx;
% keptvox_ref=keptvox;
%
%
% for eii=1:length(eis);
%     subplot(2,3,eii);
%     ei=eis(eii);
%     exp=experiments{ei};
%
%     idxs=zeros(voxn,2);
%     idxs(keptvox_ref,1)=idx_ref;
%
%     load([expdir '/' exp '/fmri/temporal/corr/' timeUnit '/voxs2voxs//LL_selfself/fc_cluster'],'C','idx','keptvox','keptT');
%     idxs(keptvox,2)=idx;
%
%     % keep only shared voxels
%     idxs(sum(idxs==0,2)~=0,:)=0;
%
%     dice=[];
%     for ki_ref=1:k;
%         for ki_test=1:k;
%             dice(ki_ref,ki_test)=2*sum(idxs(:,1)==ki_ref & idxs(:,2)==ki_test)/(sum(idxs(:,1)==ki_ref) + sum(idxs(:,2)==ki_test));
%         end
%     end
%     imagesc(dice,[0 0.5]);
%     title(exp)
% end