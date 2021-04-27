% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
% close all
loc='cluster';
set_parameters;
win=66; %(100/1.5);
timeUnit='tr';
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=25;
crop_end=20;
Fs=1/1.5;
iters=1;
network_newOrder={
    'Auditory_Language',...
    'DMN2',...
    'Attention',...
    'Executive',...
    'DMN1',...
    'Visual'...
    };
cols=jet(7);

eis=[1 2 4 11 12 13 9 10];
for eii=1%:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_HG_L.mat' ],'gdata');
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    subjs_g1=[];
    subjs_g2=[];
    for iter=1:iters;
        rng(iter)
        subjs_shuffled=randperm(listenerN);
        subjs_shuffled(ismember(subjs_shuffled,subjects_excluded{ei}))=[];
        subjs_g1(:,iter)=subjs_shuffled(1:round(length(subjs_shuffled)/2));
        subjs_g2(:,iter)=subjs_shuffled((1+round(length(subjs_shuffled)/2)):end);
    end
    
    [temp,~]=pwelch(squeeze(nanmean(gdata(:,:,1),1)),win,floor(win/2),[],Fs);
    freqN=length(temp);
    Sxx=nan(length(rnames),length(rnames),freqN,iters);
    Syy=Sxx;
    Sxy=Sxx;
    
    for sdi=1:length(rnames);
        seed=rnames{sdi};
        
        if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_' seed '.mat' ]);
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_'  seed '.mat' ],'gdata');
            gdata_seed=gdata;
            gdata_seed(:,:,subjects_excluded{ei})=NaN;
            gdata_seed=nanmean(gdata_seed,1);
            gdata_seed=zscore(gdata_seed(:,keptT,:),0,2);
            
            for tgi=1:length(rnames);
                target=rnames{tgi};
                if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_' target '.mat' ]);
                    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_'  target '.mat' ],'gdata');
                    gdata_target=gdata;
                    gdata_target(:,:,subjects_excluded{ei})=NaN;
                    gdata_target=nanmean(gdata_target,1);
                    gdata_target=zscore(gdata_target(:,keptT,:),0,2);
                    
                    for iter=1:iters;
                        x=zscore(nanmean(gdata_seed(:,:,subjs_g1(:,iter)),3),0,2);
                        y=zscore(nanmean(gdata_target(:,:,subjs_g2(:,iter)),3),0,2);
                        
                        if ~isnan(x(1));
                            [Sxy(sdi,tgi,:,iter),freq]=cpsd(x,y,win,floor(win/2),[],Fs);   % estimate Sxy
                            [Sxx(sdi,tgi,:,iter),freq]=pwelch(x,win,floor(win/2),[],Fs);   % estimate Sxx
                            [Syy(sdi,tgi,:,iter),freq]=pwelch(y,win,floor(win/2),[],Fs);    % estimate Syy
                        end
                    end
                end
            end
        end
    end
    mkdir([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/LL_gg/']);
    save([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/LL_gg/rois2rois_cpsd'  ],'rnames','keptT','freq','win','Sxx','Syy','Sxy');
end

eis=[1 2 4 11 12 13 9 10];
fsize=[45 7];
fsize=[30 27];
figure('unit','centimeter','position',[0 0 fsize]);

for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/LL_gg/rois2rois_cpsd'  ],'rnames','keptT','freq','win','Sxx','Syy','Sxy');
    subplot(3,3,eii);
    coherences=abs(Sxy).^2./(Sxx.*Syy);
    for ri=1:length(rnames);
        
        rname=rnames{ri};
        network=roi_table.network(ri);
        ni=find(ismember(network_newOrder,network));
        if ismember(ri,[1 11]);
            ciplot_claire((squeeze(coherences(ri,ri,2:end,:))'),log(freq(2:end)),cols(ni+1,:),0.3);
            hold on
        end
    end
    set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
    set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
    xlim(log([freq(2) freq(end)]))
    ylim([0 0.7])
    ylabel('Coherence');
    xlabel('Frequency (Hz)');
    %  legend(rnames_selected);
    set(gca,'fontsize',14)
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','k');
    grid on
    hold off
end
