clear all
set_parameters;
exp='sherlock';
load([ expdir exp '\sound\' exp '_listener_audenv.mat'],'aud');
load([expdir exp '/sound/onsets.mat'],'onsets_msec_wd','onsets_msec_sn','onsets_msec_pr','offsets_msec_wd','offsets_msec_sn','offsets_msec_pr','onsetsVectorMsec_wd','onsetsVectorMsec_sn','onsetsVectorMsec_pr')

pause_wd=[];
win=-5:5;
for wi=1:length(onsets_msec_wd);
     onset=round(onsets_msec_wd(wi)/1500);
  
    if onset>abs(min(win)) & (onset+max(win))<length(aud);
        pause_wd(end+1,:)=aud((onset+min(win)):(onset+max(win)));
    end;
end
pause_sn=[];
for wi=1:length(onsets_msec_sn);
       onset=round(onsets_msec_sn(wi)/1500);
  
    if onset>abs(min(win)) & (onset+max(win))<length(aud);
        pause_sn(end+1,:)=aud((onset+min(win)):(onset+max(win)));
    end;
end
pause_pr=[];
for wi=1:length(onsets_msec_pr);
     onset=round(onsets_msec_pr(wi)/1500);
   
    if onset>abs(min(win)) & (onset+max(win))<length(aud);
        pause_pr(end+1,:)=aud((onset+min(win)):(onset+max(win)));
    end;
end

cols=jet(7);
cols=cols([1 3 4 5 6 7],:)*0.95;
cols(4,:)=[1 0.85 0];

fsize=[7 8];
figure('unit','centimeter','position',[0 0 fsize]);
ciplot_claire(pause_pr,win*1.5,cols(6,:),0.3)
hold on
ciplot_claire(pause_sn,win*1.5,cols(3,:),0.3)
hold on
ciplot_claire(pause_wd,win*1.5,cols(1,:),0.3)
hold off
ylabel('Audio Envelope')
xlabel(['Time around onset (sec)'])
title([upper(exp(1)) exp(2:end)]);
set(gca,'fontsize',14)
grid on

xlim([-5 5])
ylim([min(mean(pause_pr))-80 max(mean(pause_pr))+80])
set(gca,'xtick',[-5:5])
set(gca,'xtick',[-5:5],'ytick',[]);

dur_tr=[];
pause_pr_z=zscore(pause_pr,0,2);
for pri=1:size(pause_pr,1);
temp=find(diff(pause_pr_z(pri,:)<0)~=0);
dur_tr(pri)=min(temp((temp-find(win==0)>=0)))-max(temp((temp-find(win==0)<=0)));
end
dur_masec=dur_tr*1500;
