clear all
close all
set_parameters;

ei=3;
K=28;
rname='dPCC';;
exp=experiments{ei};
load([expdir '/scripts_speaker-listener/' exp '_' rname '_hmm_findSpeakerEventInListeners.mat'],'segments_S','segments_L');
segments_S=segments_S(K,:);
segments_L=segments_L(K,:,:);
tn=size(segments_S,1);
listenerN=size(segments_L,3);

mat_l=zeros(tn,tn,listenerN);
mat=nan(tn,tn);

for s=2%:listenerN;
    for i=1:max(segments_S);
             mat_l(segments_L(:,:,s)==i,segments_L(:,:,s)==i,s)=1;
   
   %       mat_boundary=(bwboundaries(mat_l(:,:,s)));
    %      visboundaries(mat_boundary,'color','b');
    end
end


mat_s=zeros(tn,tn);
for i=1:max(segments_S);
    mat(segments_S==i,segments_S==i)=i;
    mat_s(segments_S==i,segments_S==i)=1;
end
mat_boundary=(bwboundaries(mat_s));


figure;
hold on
imagesc(mean(mat_l,3)); colormap gray;
visboundaries(mat_boundary,'color','r');
hold off

% figure
%     plot([1 K],[1 K],'r')
%     hold on
%     
% for s=1:listenerN;
%      scatter(squeeze(segments_S),squeeze(segments_L(:,:,s)),20,'filled','k','MarkerFaceAlpha',1/18);
% 
% end
% hold off
% 
% 
% 
%      
%      
     
    

