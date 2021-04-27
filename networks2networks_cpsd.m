% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
% nonzero coherence at frequency is becasue  of the welch method and could
% vry with window length
clear all
% close all
% loc='cluster';
set_parameters;

win=100; %(100/1.5);
eis=[1 2 4 11 12 13 9 10];
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
crop_start=25;
crop_end=20;
win=80;

for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/vox2vox/LL_leave1out/isc'   ],'peakLags','peaks','keptvox');
    maskPeakLag0=zeros(voxn,1);
    thr=sort(peaks(peakLags==0),'descend');
    thr=thr(round(length(keptvox)*0.3));
    maskPeakLag0(keptvox(peaks>thr & peakLags'==0))=1;
    clear peakLags peaks keptvox
    
    Fs=1/tr(ei);
    
    for sdi=1:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox');
        % keep only voxel with peakLag 0
        disp(sum(ismember(keptvox,find(maskPeakLag0))))
        gdata_seed=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
        gdata_seed(:,:,subjects_excluded{ei})=NaN;
        gdata_seed=nanmean(gdata_seed,1);
        
        [~,tn,listenerN]=size(gdata_seed);
        keptT=(crop_start+1):(tn-crop_end);
        lags=-floor((length(keptT)-1)/2):floor((length(keptT)-1)/2);
        
        r=nan([length(networks)  length(lags) listenerN ]);
        
        for tgi=1:size(networks);
            network=networks{tgi};
            
            if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ]);
                load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' network '.mat' ],'gdata','keptvox');
                % keep only voxels with peakLag 0
                gdata=gdata(ismember(keptvox,find(maskPeakLag0)),:,:);
                gdata(:,:,subjects_excluded{ei})=NaN;
                gdata=nanmean(gdata,1);
                
                for si=1:listenerN;
                    othersi=1:listenerN;
                    othersi=othersi(othersi~=si);
                    
                    y=zscore(nanmean(zscore(gdata(:,keptT,othersi),0,2),3),0,2);
                    x=zscore(gdata_seed(:,keptT,si),0,2);
                    
                    [Sxy(sdi,tgi,:,si),freq]=cpsd(x,y,win,win/2,[],Fs);   % estimate Sxy
                    [Sxx,freq]=pwelch(x,win,win/2,[],Fs);   % estimate Sxx
                    [Syy,freq]=pwelch(y,win,win/2,[],Fs);    % estimate Syy
                    
                end
            end
        end
    end
    save([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_cpsd'  ],'Sxy','networks','keptT','freq');
end



fsize=[30 18];
figure('unit','centimeter','position',[0 0 fsize],'color','k');
nis=[2     4     1     5     3     6];
cols=jet(7);
eis=[9 10];
for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_cpsd'  ],'Sxy','networks','keptT','freq');
    subplot(2,3,eii);
    
    for tgi=1:6;
        ciplot_claire(squeeze(abs(Sxy(nis(tgi),nis(tgi),:,:)))',freq,cols(tgi+1,:),0);
        
        hold on
    end
    xlim([0 0.1]);
    
    % ylim([-5 25]);
    gray=[0.7 0.7 0.7];
    set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
    
    line([0 0 ],get(gca,'ylim'),'color',gray);
    line(get(gca,'xlim'),[0 0],'color',gray);
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','w');
    
    grid on
    hold off
end


%%
fsize=[30 18];
figure('unit','centimeter','position',[0 0 fsize],'color','k');
nis=[2     4     1     5     3     6];
cols=jet(7);
for eii=7:8;%1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/network2networks/' froidir '/LL_leave1out/networks2networks_cpsd'  ],'Sxy','networks','keptT','freq');
    subplot(2,3,eii);
    
    mx=abs(nanmean(Sxy(2,:,:,:),4));
    mx=max(mx(:));
    for tgi=1:6;
        %   ciplot_claire(squeeze(abs(Sxy(3,nis(tgi),:,:)))',freq,cols(tgi+1,:),0);
        
        temp=squeeze(nanmean(Sxy(nis(tgi),nis(tgi),:,:),4));
        scatter(nanmean(real(temp),2)/mx,nanmean(imag(temp),2)/mx,15,cols(tgi+1,:),'filled');
        % % scatter(cos(angle(temp)),sin(angle(temp)),40,cols(tgi+1,:),'markerfacealpha',0.4);
        hold on
    end
    xlim([-1.3 1.3]);
    ylim([-1.3 1.3]);
    gray=[0.7 0.7 0.7];
    set(gca,'color','k','XColor',gray,'YColor',gray,'GridColor',gray)
    viscircles([0 0],1,'color',gray,'linewidth',0.2);
    line([0 0 ],get(gca,'ylim'),'color',gray);
    line(get(gca,'xlim'),[0 0],'color',gray);
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','w');
    % xlim([min(freq) max(freq)])
    % xlim([0 0.1])
    % line([0.04 0.04],get(gca,'ylim'),'color','k')
    
    hold off
end
