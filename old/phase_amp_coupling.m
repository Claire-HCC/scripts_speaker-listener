sdi=2;
tgi=3;
seed=networks{sdi};
load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox');
% keep only voxel with peakLag 0
% disp(sum(ismember(keptvox,find(maskPeakLag0))))
gdata_seed=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
gdata_seed(:,:,subjects_excluded{ei})=NaN;
gdata_seed=nanmean(gdata_seed,1);
%    gdata_seed=gdata_seed(:,keptT,:);
%   gdata_seed=gdata_seed(:,:,subjs_g1);
%  x=nanmean(zscore(nanmean(gdata_seed,1),0,2),3)';
lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
target=networks{tgi};
load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' target '.mat' ],'gdata','keptvox');
% keep only voxel with peakLag 0
% disp(sum(ismember(keptvox,find(maskPeakLag0))))
[~,tn,listenerN]=size(gdata);
keptT=(crop_start+1):(tn-crop_end);
gdata_target=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
gdata_target(:,:,subjects_excluded{ei})=NaN;
gdata_target=nanmean(gdata_target,1);
y=nanmean(zscore(gdata_target(:,keptT,:),0,2),3);
x=nanmean(zscore(gdata_seed(:,keptT,:),0,2),3);


fi=2;
xf = filter1(filtertypes{fi},double(x),'fc',cutoff,'fs',Fs);
yf = filter1(filtertypes{find(~ismember([1:2],fi))},double(y),'fc',cutoff,'fs',Fs);
plot([xf; yf]')
xf_env=abs(hilbert(xf));

plot(lags,lagcorr(xf_env ,yf,lags))
scatter(angle(hilbert(yf)),xf_env)
bar(rad2deg(-3.3:0.33:3.3),grpstats(xf_env,yf_pha_bined))