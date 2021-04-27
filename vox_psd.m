% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
close all

loc='cluster';
set_parameters;
win=66; %(100/1.5);

timeUnit='tr';
crop_start=25;
crop_end=20;

eis=[1 2 4 11 12 9 10 13];
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain/listenerAll_iscmasked.mat' ],'gdata','keptvox');
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    [temp,~]=pwelch(squeeze(nanmean(gdata(:,:,1),1)),win,floor(win/2),[],1/tr(ei));
    freqN=length(temp);
    Sxx=nan(length(keptvox),freqN);
    
    gdata(:,:,subjects_excluded{ei})=NaN;
    gdata=gdata(:,keptT,:);
    % so that subject with larger mean signals would not gain larger
    % weighting
    gdata=gdata-nanmean(gdata,2);
    % average across subjects first. In Stephens (2013), the results were more salient with this step, presumablly due to the attenuation of individual subject's intrnsic signals.
    gdata=nanmean(gdata,3);
    % Since the variance is equal across all voxels, all spectra have the same integrated area.
    gdata=zscore(gdata,0,2);
    
    tic
    for vi=1:length(keptvox)
        [Sxx(vi,:) freq] = pwelch(gdata(vi,:)',win,win/2,[],1/tr(ei));
    end
    toc 
     mkdir([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/wholeBrain/L_g/']);
    save([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/wholeBrain/L_g/wholeBrain_psd'  ],'Sxx','keptvox','keptT','freq','win');
end
