clear all
set_parameters;
exp='sherlock';
load([ expdir exp '\sound\' exp '_listener_audenv_msec.mat'],'aud');
load([expdir exp '/sound/onsets.mat'],'onsets_msec_wd','onsets_msec_sn','onsets_msec_pr','offsets_msec_wd','offsets_msec_sn','offsets_msec_pr','onsetsVectorMsec_wd','onsetsVectorMsec_sn','onsetsVectorMsec_pr')

% smooth the audio=evelope. Otherwise it ould be difficult to decide the
% length of between-paragraph boundary.
aud=smooth(aud',1000);

pause_wd=[];
win=-6000:6000;
for wi=1:length(onsets_msec_wd);
    %  onset=round(onsets_msec_wd(wi)/1500);
    onset=onsets_msec_wd(wi);
    if onset>abs(min(win)) & (onset+max(win))<length(aud);
        pause_wd(end+1,:)=aud((onset+min(win)):(onset+max(win)));
    end;
end
pause_sn=[];
for wi=1:length(onsets_msec_sn);
    %   onset=round(onsets_msec_sn(wi)/1500);
    onset=onsets_msec_sn(wi);
    if onset>abs(min(win)) & (onset+max(win))<length(aud);
        pause_sn(end+1,:)=aud((onset+min(win)):(onset+max(win)));
    end;
end
pause_pr=[];
for wi=1:length(onsets_msec_pr);
    %  onset=round(onsets_msec_pr(wi)/1500);
    onset=onsets_msec_pr(wi);
    if onset>abs(min(win)) & (onset+max(win))<length(aud);
        pause_pr(end+1,:)=aud((onset+min(win)):(onset+max(win)));
    end;
end

fsize=[12 12];
figure('unit','centimeter','position',[0 0 fsize]);
ciplot_claire(pause_wd,win/1000,'b',0.3)
hold on
ciplot_claire(pause_sn,win/1000,'g',0.3)
hold on
ciplot_claire(pause_pr,win/1000,'r',0.3)
hold off
ylabel('Audio Envelope')
xlabel(['Time around onset (sec)'])
title(exp)
set(gca,'fontsize',14)
grid on

xlim([min(win) max(win)]/1000)
ylim([920 1080])
set(gca,'xtick',[-5:5],'ytick',[])


dur_boundary=[];
pause_pr_z=zscore(pause_pr,0,2);
for pri=1:size(pause_pr,1);
    temp=find(diff(pause_pr_z(pri,:)<0)~=0);
    dur_boundary(pri)=min(temp((temp-find(win==0)>=0)))-max(temp((temp-find(win==0)<=0)));
end
save([expdir exp '/sound/paragraphBoundaryDuration_msec.mat'],'dur_boundary');
fsize=[9 9];
figure('unit','centimeter','position',[0 0 fsize]);
histogram(dur_boundary/1000,10,'facecolor','k','binLimits',[0 max(win)/1000]);
xlabel('Boundary duration (sec)');
ylabel('Paragraph N')
title([upper(exp(1)) exp(2:end)])
set(gca,'fontsize',12);
xlim([0 max(win)]/1000)

fsize=[10 8];
figure('unit','centimeter','position',[0 0 fsize]);
imagesc(pause_pr)
colormap hot;
set(gca,'xtick',1:1000:length(win),'xticklabel',[(min(win)/1000):(max(win)/1000)],'ytick','');
xlabel('Time from paragraph onset (sec)')
ylabel('Paragraph')
line([find(win==0) find(win==0)],[-0.5 21.5],'color','w','linewidth',3)
title([upper(exp(1)) exp(2:end)])
set(gca,'fontsize',12)
c=colorbar;
c.Label.String='Audio-envelope'
c.Ticks=[];