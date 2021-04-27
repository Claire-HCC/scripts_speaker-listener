% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
clear all
% close all
% loc='cluster';
set_parameters;
win=66; %(100/1.5);
timeUnit='tr';
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
crop_start=25;
crop_end=20;
Fs=1/1.5;
network_newOrder={
    'Auditory_Language',...
    'DMN2',...
    'Attention',...
    'Executive',...
    'DMN1',...
    'Visual'...
    };
cols=jet(7);
%
% eis=[1 2 4 11 12 13 9 10];
% for eii=1:length(eis);
%     ei=eis(eii);
%     exp=experiments{ei};
%
%     load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_HG_L.mat' ],'gdata');
%     [~,tn,listenerN]=size(gdata);
%     keptT=(crop_start+1):(tn-crop_end);
%
%     [temp,~]=pwelch(squeeze(nanmean(gdata(:,:,1),1)),win,floor(win/2),[],Fs);
%     freqN=length(temp);
%     Sxx=nan(length(rnames),length(rnames),freqN,listenerN);
%     Syy=Sxx;
%     Sxy=Sxx;
%
%     for sdi=1:length(rnames);
%         seed=rnames{sdi};
%         if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_' seed '.mat' ]);
%             load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_'  seed '.mat' ],'gdata');
%             gdata_seed=gdata;
%             gdata_seed(:,:,subjects_excluded{ei})=NaN;
%             gdata_seed=nanmean(gdata_seed,1);
%             gdata_seed=zscore(gdata_seed(:,keptT,:),0,2);
%
%             for tgi=1:length(rnames);
%                 target=rnames{tgi};
%                 if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_' target '.mat' ]);
%                     load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_iscmasked_'  target '.mat' ],'gdata');
%                     gdata_target=gdata;
%                     gdata_target(:,:,subjects_excluded{ei})=NaN;
%                     gdata_target=nanmean(gdata_target,1);
%                     gdata_target=zscore(gdata_target(:,keptT,:),0,2);
%
%                     for si=1:listenerN;
%                         othersi=1:listenerN;
%                         othersi=othersi(othersi~=si);
%
%                         y=zscore(nanmean(gdata_target(:,:,othersi),3),0,2);
%                         x=zscore(gdata_seed(:,:,si),0,2);
%
%                         if ~isnan(x(1));
%                             [Sxy(sdi,tgi,:,si),freq]=cpsd(x,y,win,floor(win/2),[],Fs);   % estimate Sxy
%                             [Sxx(sdi,tgi,:,si),freq]=pwelch(x,win,floor(win/2),[],Fs);   % estimate Sxx
%                             [Syy(sdi,tgi,:,si),freq]=pwelch(y,win,floor(win/2),[],Fs);    % estimate Syy
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     delete([expdir '/' exp '/fmri/temporal/circularlagcorr/' timeUnit '/roi2rois/' froidir '/LL_leave1out/rois2rois_coherence'  ]);
%     save([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/LL_leave1out/rois2rois_cpsd'  ],'rnames','keptT','freq','win','Sxx','Syy','Sxy');
% end

eis=[1 2 4 11 12 13 9 10];

fsize=[30 27];
figure('unit','centimeter','position',[0 0 fsize]);

for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/roi/' froidir '/LL_leave1out/rois2rois_cpsd'  ],'rnames','keptT','freq','win','Sxx','Syy','Sxy');
    
    subplot(3,3,eii);
    coherences=abs(Sxy).^2./(Sxx.*Syy);
    temp=squeeze(coherences(11,11,:,:))'-squeeze(coherences(1,1,:,:))';
    ciplot_claire(temp(:,2:end),log(freq(2:end)),'k',0.3);
    
    %     for ri=1:length(rnames);
    %         rname=rnames{ri};
    %         network=roi_table.network(ri);
    %         ni=find(ismember(network_newOrder,network));
    %         if ismember(ri,[1 11]);
    %            ciplot_claire((squeeze(coherences(nis(ni),nis(ni),2:end,:))'),log(freq(2:end)),cols(ni+1,:),0.3);
    %             %    plot(log(freq(2:end)),squeeze(nanmean(coherences(ri,ri,2:end,:),4)),'color',cols(ni+1,:));
    %             hold on
    %         end
    %     end
    
    set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
    set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
    xlim(log([freq(2) freq(end)]))
    %   ylim([0 40])
    %  ylabel('Coherence');
    xlabel('Frequency (Hz)');
    %  legend(rnames_selected);
    set(gca,'fontsize',14)
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')],'color','k');
    %¡@grid on
    ylabel({'Coheherence','vPCUN > HG L'})
    line([get(gca,'xlim')],[0 0 ],'color','k','linewidth',0.5)
end
