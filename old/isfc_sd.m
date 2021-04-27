function isfc_sd(ei)

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;
iters=1;

exp=experiments{ei};
mkdir([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/voxs2voxs/LL_leave1out/']);

load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/vox2vox//LL_leave1out/lag-20-20'  ],'peakLags','peaks','keptvox');
peakLag0_mask=zeros(voxn,1);
% keep vox with top 30 isc, within reginos with peak lag 0
thr=sort(peaks(peakLags==0),'descend');
thr=thr(round(length(keptvox)*0.3));
peakLag0_mask(keptvox(peaks>thr & peakLags==0))=1;

load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll.mat' ],'gdata','keptvox');

gdata=gdata(ismember(keptvox,find(peakLag0_mask)),:,:);
keptvox=keptvox(ismember(keptvox,find(peakLag0_mask)));

gdata(:,:,subjects_excluded{ei})=NaN;
[~,tn,listenerN]=size(gdata);

keptT=(crop_start+1):(tn-crop_end);
lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
subjs_g1=[];
subjs_g2=[];
for iter=1:iters;
    rng(iter)
    subjs_shuffled=randperm(listenerN);
    subjs_shuffled(ismember(subjs_shuffled,subjects_excluded{ei}))=[];
    subjs_g1(:,iter)=subjs_shuffled(1:round(length(subjs_shuffled)/2));
    subjs_g2(:,iter)=subjs_shuffled((1+round(length(subjs_shuffled)/2)):end);
end


g1=nanmean(zscore(gdata(:,keptT,subjs_g1(:,1)),0,2),3);
g2=nanmean(zscore(gdata(:,keptT,subjs_g2(:,1)),0,2),3);

sd=nan(length(keptvox),length(keptvox));
peaks=nan(length(keptvox),length(keptvox));
peakLags=nan(length(keptvox),length(keptvox));
for vi=1:length(keptvox);
    disp(vi1)
    
    temp=(circularlagcorr_byCol(repmat(g1(vi,:)',1,length(keptvox)),g2',lags));
    sd(vi,:)=std(temp);
    [~,lagi]=max(abs(temp),[],1);
    lagind=sub2ind(size(temp),lagi,(1:size(temp,2)));
    peaks(vi,:)=temp(lagind);
    peakLags(vi,:)=lags(lagi);
end

save([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/voxs2voxs/LL_leave1out/isfc' '_circularlag' num2str(min(lags)) '-' num2str(max(lags)) '_peaks.mat' ],'peaks','peakLags','sd','keptvox','subjs_g1','subjs_g2','-v7.3');
