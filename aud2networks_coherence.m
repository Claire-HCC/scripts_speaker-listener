% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
% nonzero coherence at frequency is becasue  of the welch method and could
% vry with window length
clear all
% close all
% loc='cluster';
set_parameters;

win=100; %(100/1.5);
eis=[1 2 4 11 12];
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;
win=100;

nis=[2     4     1     5     3     6];
cols=jet(7);
for ei=[1 2 4 11 12]
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    Fs=1/tr(ei);
    
    load([ expdir exp '\sound\' exp '_listener_audhrf.mat'],'aud');
    
    for tgi=1:length(networks);
        target=networks{tgi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' target '.mat' ],'gdata','keptvox');
        % keep only voxel with peakLag 0
        disp(sum(ismember(keptvox,find(maskPeakLag0))))
        gdata_target=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
        gdata_target(:,:,subjects_excluded{ei})=NaN;
        gdata_target=nanmean(gdata_target,1);
        
        [~,tn,listenerN]=size(gdata_target);
        keptT=(crop_start+1):(tn-crop_end);
        lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
        
        y=zscore(nanmean(zscore(gdata_target(:,keptT,:),0,2),3),0,2);
        x=zscore(aud(keptT));
        
        [Sxy,freq]=cpsd(x,y,win,win/2,[],Fs);   % estimate Sxy
        [Sxx,freq]=pwelch(x,win,win/2,[],Fs);   % estimate Sxx
        [Syy,freq]=pwelch(y,win,win/2,[],Fs);    % estimate Syy
        % [Cms,freq]=mscohere(   x,y,win,win/2,[],Fs);    % estimate ordinary
        Cxy = Sxy ./ sqrt( Sxx .* Syy );
        coh(tgi,:)=abs(Cxy).^2;
        Kxy  = real( Sxy );
        Qxy  = imag( Sxy );
        pha(tgi,:)  = atan2( Qxy, Kxy );
    end
     save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/aud2networks_coherence'  ],'pha','coh','networks','keptT','freq');
 
end



nis=[2     4     1     5     3     6];
cols=jet(7);

for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/aud2networks_coherence'  ],'pha','coh','networks','keptT','freq');
    subplot(2,3,eii);
    
    for tgi=1:6;
        plot(freq,pha(nis(tgi),:)','color',cols(tgi+1,:),'linewidth',2);
        hold on
    end
    title([strrep(exp,'_',' ')]);
    xlim([min(freq) max(freq)])
    %   xlim([0 0.1])
    grid on
    hold off
end

