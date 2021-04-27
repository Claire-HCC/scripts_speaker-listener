
clear all
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;
cutoff=[0.04];
Fs=1/1.5;
filtertypes={'lp','hp'};

for ei=[13];%[1 4 11 12 9 10];
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seeds{1} '.mat' ],'gdata','keptvox');
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
    r=nan(length(networks),length(networks),length(lags),length(filtertypes),1);
    
    gdataf=[];
    for sdi=1:length(networks);
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
        
        for tgi=1:6;
            target=networks{tgi};
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' target '.mat' ],'gdata','keptvox');
            % keep only voxel with peakLag 0
            % disp(sum(ismember(keptvox,find(maskPeakLag0))))
            [~,tn,listenerN]=size(gdata);
            keptT=(crop_start+1):(tn-crop_end);
            
            gdata_target=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
            gdata_target(:,:,subjects_excluded{ei})=NaN;
            gdata_target=nanmean(gdata_target,1);
            
            for fi=1:2;
                filtertype=filtertypes{fi};
                
                gdata_targetf=[];
                gdata_seedf=[];
                for si=1:listenerN;
                    gdata_targetf(:,:,si)=filter1(filtertype,double(gdata_target(:,:,si)),'fc',cutoff,'fs',Fs);
                    gdata_seedf(:,:,si)=filter1(filtertype,double(gdata_seed(:,:,si)),'fc',cutoff,'fs',Fs);
                end
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    y=nanmean(zscore(gdata_targetf(:,keptT,othersi),0,2),3);
                    x=zscore(gdata_seedf(:,keptT,si),0,2);
                    
                    r(sdi,tgi,:,fi,si)=circularlagcorr(x,y,lags)';
                    
                end
            end
            gdataf(sdi,:,:,fi)= gdata_targetf;
        end
        
    end
    
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_filtered' '_cutoff'  strrep(num2str(cutoff),'0.','')  ],'cutoff','r','networks','keptT','filtertypes','gdataf');
end


nis=[2     4     1     5     3     6];
cols=jet(7);
eis=[1 2 4 11 12 ];
for fi=1:2;
    fsize=[35 8];
    figure('unit','centimeter','position',[0 0 fsize]);%,'color','k');
    for eii=1%:length(eis);
        ei=eis(eii);
        exp=experiments{ei};
        % load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_filtered'  ],'r','networks','keptT','filtertypes');
        rzm=nanmean(atanh(r),5);
        lags=(size(rzm,3)-1)/2;lags=[-lags:lags];
        subplot(1,5,eii)
        for tgi=1:6;
            plot(lags,squeeze(rzm(2,nis(tgi),:,fi)),'linewidth',2,'color',cols(tgi+1,:)); hold on;
            xlim([-20 20]); ylim([-0.3 0.5]);
            grid on
        end
        title([upper(exp(1)) strrep(exp(2:end),'_',' ')])
    end
end

sdi=3;
for eii=1%:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    figure;
    % load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_filtered'  ],'r','networks','keptT','filtertypes');
    rzm=nanmean(atanh(r),5);
    lags=(size(rzm,3)-1)/2;lags=[-lags:lags];
    for fi=1:2;
        subplot(1,2,1)
        imagesc(squeeze(rzm(sdi,nis,:,1)),[-0.2 0.2]);
        colormap jet
        subplot(1,2,2)
        imagesc(squeeze(rzm(sdi,nis,:,2)),[-0.2 0.2]);
        colormap jet
    end
    title(exp)
end


