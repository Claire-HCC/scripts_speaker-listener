clear all

loc='cluster';
set_parameters;
timeUnit='tr' ;
crop_start=25;
crop_end=20;
froidir='isc_peak';
rois={'HG_L','Angular_L','precuneus'};
perc=0.30;

for ei=1:10;
    exp=exp_parameters.experiments{ei};
    mkdir([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_leave1out/']);
    
    f= sprintf('%s/%s/fmri/timeseries/%s/wholeBrain/listenerAll_isc%dPercMasked.mat',expdir,exp,timeUnit,perc*100);
    load(f,'gdata','keptvox');
    gdata_target=gdata;
    gdata_target(:,:,exp_parameters.subjects_excluded{ei})=NaN;
    
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    r=nan(length(keptvox),length(lags),listenerN);
    
    for ri=1;%1:length(rois);
        roi=rois{ri};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' roi '.mat' ],'gdata');
        gdata_seed= nanmean(gdata,1);
        gdata_seed(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        
        for si=1:listenerN;
            othersi=1:listenerN;
            othersi=othersi(~ismember(othersi,si));
            
            y=nanmean(zscore(gdata_target(:,keptT,othersi),0,2),3)';
            x=nanmean(zscore(gdata_seed(:,keptT,si),0,2),3)';
            x=repmat(x,1,size(y,2));
            
            r(:,:,si)=circularlagcorr_byCol(x,y,lags)';
        end
        
        save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2vox/' froidir '/LL_leave1out/' roi '.mat'   ],'r','lags','keptvox','keptT','-v7.3');
    end
end
