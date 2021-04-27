% find the peak nearest to lag 0 instead of the absolute peak
clear all

loc='mypc';
set_parameters;

timeUnit='tr' ;
froidir='isc_peak';
lags=-15:15;
perc=0.30;

rzms=[];
rzms_timeReversed=[];
rois={'HG_L','Angular_L','precuneus'};

exp='crossExps';
for ri=1;%1:length(rois);
    roi=rois{ri};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir  '/LL_selfself/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReverse_fdr_peakLags.mat' ],...
        '','rzms','rzms_timeReversed','rzmsm','rzmsm_timeReversed','lags','peakLags','peaks','rzm','pfdr_peaks','pfdr_exps_peaks','pfdr_time_peaks','tmid','npeaks');
    
    sig_fc=(pfdr_peaks<.05 & (abs(npeaks)<peaks | isnan(npeaks) ));
    peakLags_fc=peakLags;
    
     load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir  '/LL_leave1out/' roi '2voxIsc' num2str(perc*100) 'PercMasked_lag' num2str(min(lags)) '-' num2str(max(lags)) '_timeReverse_fdr_peakLags.mat' ],...
        '','rzms','rzms_timeReversed','rzmsm','rzmsm_timeReversed','lags','peakLags','peaks','rzm','pfdr_peaks','pfdr_exps_peaks','pfdr_time_peaks','tmid','npeaks');
    
    sig_isfc=(pfdr_peaks<.05 & (abs(npeaks)<peaks | isnan(npeaks) ));
    peakLags_isfc=peakLags;
    
   
end
sig=(sig_fc & sig_isfc);

x=peakLags_fc(sig)*1.5;
y=peakLags_isfc(sig)*1.5;

t = 2*pi*rand(length(x),1);
r = 0.75*sqrt(rand(length(x),1));
 
            x_jittered = x + r.*cos(t);
            y_jittered = y + r.*sin(t);
            scatter(x_jittered, y_jittered,5,[0.5 0.5 0.5 ],'filed','markerfacealpha',0.3)
            xlim([min(y)-1 max(y)+1]);
              ylim([min(y)-1 max(y)+1])