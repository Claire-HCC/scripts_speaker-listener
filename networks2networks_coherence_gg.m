% the effect of window length
% https://sapienlabs.org/factors-that-impact-power-spectrum-density-estimation/
% nonzero coherence at frequency is becasue  of the welch method and could
% vry with window length
clear all
close all
% loc='cluster';
set_parameters;

eis=[ 1 2  4 9:13  ];
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
networks=unique(table2array(roi_table(:,2)));
seeds=networks;
crop_start=25;
crop_end=20;
win=66; %(100/1.5);
Fs=1/1.5;

iters=100;;
% for eii=1:length(eis);
%     ei=eis(eii);
%     exp=experiments{ei};
%     
%     load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_iscmasked_' networks{1} '.mat' ],'gdata','keptvox');
%     [~,tn,listenerN]=size(gdata);
%     keptT=(crop_start+1):(tn-crop_end);
%     
%     subjs_g1=[];
%     subjs_g2=[];
%     for iter=1:iters;
%         rng(iter)
%         subjs_shuffled=randperm(listenerN);
%         subjs_shuffled(ismember(subjs_shuffled,subjects_excluded{ei}))=[];
%         subjs_g1(:,iter)=subjs_shuffled(1:round(length(subjs_shuffled)/2));
%         subjs_g2(:,iter)=subjs_shuffled((1+round(length(subjs_shuffled)/2)):end);
%     end
%     
%     [temp,~]=pwelch(squeeze(nanmean(gdata(:,:,1),1)),win,floor(win/2),[],Fs);
%     freqN=length(temp);
%     Sxx=nan(length(networks),length(networks),freqN,listenerN);
%     Syy=Sxx;
%     Sxy=Sxx;
%     
%     for sdi=1:length(networks);
%         seed=networks{sdi};
%         load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_iscmasked_' seed '.mat' ],'gdata','keptvox');
%         gdata_seed=gdata;
%         gdata_seed(:,:,subjects_excluded{ei})=NaN;
%         gdata_seed=nanmean(gdata_seed,1);
%         gdata_seed=zscore(gdata_seed(:,keptT,:),0,2);
%         
%         for tgi=1:size(networks);
%             network=networks{tgi};
%             
%             if exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_iscmasked_' network '.mat' ]);
%                 load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_iscmasked_' network '.mat' ],'gdata','keptvox');
%                 gdata_target=gdata;
%                 gdata_target(:,:,subjects_excluded{ei})=NaN;
%                 gdata_target=nanmean(gdata_target,1);
%                 gdata_target=zscore(gdata_target(:,keptT,:),0,2);
%                 
%                 for iter=1:iters;
%                     
%                     y=zscore(nanmean(zscore(gdata_target(:,:,subjs_g1(:,iter)),0,2),3),0,2);
%                     x=zscore(nanmean(zscore(gdata_seed(:,:,subjs_g2(:,iter)),0,2),3),0,2);
%                     
%                     [Sxy(sdi,tgi,:,iter),freq]=cpsd(x,y,win,floor(win/2),[],Fs);   % estimate Sxy
%                     [Sxx(sdi,tgi,:,iter),freq]=pwelch(x,win,floor(win/2),[],Fs);   % estimate Sxx
%                     [Syy(sdi,tgi,:,iter),freq]=pwelch(y,win,floor(win/2),[],Fs);    % estimate Syy
%                     
%                 end
%             end
%         end
%     end
%     mkdir([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/network/' froidir '/LL_gg/'])
%     save([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/network/' froidir '/LL_gg/networks2networks_cpsd' ],'networks','keptT','freq','win','Sxx','Syy','Sxy');
%     
% end

eis=[ 1 2  4 11 12 13 9 10  ];
fsize=[30 27];
figure('unit','centimeter','position',[0 0 fsize]);
nis=[2     4     1     5     3     6];
cols=jet(7);

for eii=1:length(eis);
    ei=eis(eii);
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal/frequency/' timeUnit '/network/' froidir '/LL_gg/networks2networks_cpsd' ],'networks','keptT','freq','win','Sxx','Syy','Sxy');
    subplot(3,3,eii);
    
    coherences=abs(Sxy).^2./(Sxx.*Syy);
    
    for ni=1:length(networks);
        ciplot_claire((squeeze(coherences(nis(ni),nis(ni),:,:))'),log(freq),cols(ni+1,:),0.1);
        hold on
    end
    
    title([upper(exp(1)) strrep(exp(2:end),'_',' ')])
    xlabel('Frequency (Hz)');
    
    %  ylim([0 50])
    set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
    set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
    xlim(log([0.01 0.33]))
    grid on
    set(gca,'fontsize',14);
    ylabel('Coherence')
    hold off
end
