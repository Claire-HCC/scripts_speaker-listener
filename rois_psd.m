% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
% close all
set_parameters;
win=66; %(100/1.5);

timeUnit='tr';
froidir='aal';
rnames={'HG_L','Angular_L','Precuneus_L'};


crop_start=25;
crop_end=20;
Fs=1/1.5;

for ei=1:10;
    exp=exp_parameters.experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/L_g/']);
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_HG_L.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    [temp,~]=pwelch(squeeze(nanmean(gdata(:,:,1),1)),win,floor(win/2),[],Fs);
    freqN=length(temp);
    Sxx=nan(length(rnames),freqN);
    
    for ri=1:length(rnames);
        
        rname=rnames{ri};
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_'  rname '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_'  rname '.mat' ],'gdata');
             %% stephens (2013) method
        gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        gdata=gdata(:,keptT,:);
        gdata=zscore(gdata,0,2);
        % average across subjects first. In Stephens (2013), the results were more salient with this step, presumablly due to the attenuation of individual subject's intrnsic signals.
        gdata=nanmean(gdata,3);
        % Stephens: Since the variance is equal across all voxels, all
        % spectra have the same integrated area. CHC:sum(Sxx_temp)*(freq(2)-freq(1));  the integrted areas
        % would only be exactly the same using periodogram: 
         gdata=zscore(gdata,0,2);
        [Sxx_temp freq] = pwelch(gdata',win,win/2,[],1/exp_parameters.tr(ei));
        
        
        Sxx(ri,:)=nanmean(Sxx_temp,2);
        end
    end

    save([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/L_g/rois_psd'  ],'Sxx','rnames','keptT','freq','win');
end
