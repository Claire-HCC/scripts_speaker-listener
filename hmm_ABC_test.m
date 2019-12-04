clear all
close all
set_parameters;

load([expdir '..\ABC\texts\segment_time.mat'])
segmentv_inTr2=zeros(size(segmentv_inTr));

ei=4;
rname='ABCshen222';;
exp=experiments{ei};
load([expdir '/scripts_speaker-listener/' exp '_' rname '_hmm.mat'],'segments');
segmentv_inTr_hmm=round(mean(segments(18,:,:),3));

segmentv_inTr_text=zeros(size(segmentv_inTr));
n=1;
for i=1:45;
    segmentv_inTr_text(segment_tr(i,1):segment_tr(i,2))=n;
    n=n+1;
    if i<45;
    segmentv_inTr_text((segment_tr(i,2)+1):((segment_tr(i+1,1)-1)))=n;
    n=n+1;
    end
end

mat_text=zeros(2225,2225);
for i=1:max(segmentv_inTr_text);
    mat_text(segmentv_inTr_text==i,segmentv_inTr_text==i)=1;
end

mat_hmm=zeros(2225,2225);
for i=1:max(segmentv_inTr_hmm);
    mat_hmm(segmentv_inTr_hmm==i,segmentv_inTr_hmm==i)=1;
end

mat_text_boundary=(bwboundaries(mat_text));


figure;
imagesc(mat_hmm);
hold on
visboundaries(mat_text_boundary);
hold off
