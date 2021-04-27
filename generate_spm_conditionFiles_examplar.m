clear all
set_parameters;

modelName='model_examplar';
exp='sherlock';

names={'word','sentence','paragraph'};

onsets=cell(1,length(names));

durations=cell(1,length(names));

sn_wdN=10;
pr_snN=5;
dur=0.3;
wdN=(sn_wdN*pr_snN);
onsets{1}=[3:dur:(3+dur*(wdN-1))];
onsets{2}=onsets{1}((0:sn_wdN:(wdN-1))+1);
onsets{3}=onsets{1}((0:wdN:(wdN-1))+1);


durations{1}=(dur-0.01)*ones(size(onsets{1}));
durations{2}=(dur*sn_wdN-0.01)*ones(size(onsets{2}));
durations{3}=(dur*wdN-0.01)*ones(size(onsets{3}));


pmod(1).name{1}='sn_wdN';
pmod(2).name{1}='pr_snN';


pmod(1).param{1} = log(repmat(1:sn_wdN,1,pr_snN)+1);
pmod(1).poly{1}  = 1;
pmod(2).param{1} = log([1:pr_snN]+1);
pmod(2).poly{1}  = 1;


mkdir([expdir exp '\fmri\temporal\stimulus2brain_glm\' modelName '\']);
save([expdir exp '\fmri\temporal\stimulus2brain_glm\' modelName '\conditionSpec.mat'],'onsets','durations','names','pmod');


figure;
tn=max(find(SPM.Sess.U(2).u(:,2)>0));
subplot(2,2,1);
plot(SPM.Sess.U(2).u(:,2),'k', 'linewidth',2);
set(gca,'xtick',[],'ytick',[],'fontsize',14);
xlim([0 tn+300]);
title('Neural activity')
ylabel({'Paragraph level','accumulatio by sentence'});
subplot(2,2,2);
title('BOLD response')
plot(SPM.xX.X(:,4),'k', 'linewidth',2);
set(gca,'xtick',[],'ytick',[],'fontsize',14);
xlim([0 round(tn/150)+10])
subplot(2,2,3);
plot(SPM.Sess.U(1).u(:,2),'k', 'linewidth',2);
set(gca,'xtick',[],'ytick',[],'fontsize',14);
xlim([0 tn+300]);
ylabel({'Sentence level','accumulatin by word'});
subplot(2,2,4);
plot(SPM.xX.X(:,2),'k', 'linewidth',2);
set(gca,'xtick',[],'ytick',[],'fontsize',14);
xlim([0 round(tn/150)+10])

 figure; 
 lags=-9:9;
 r=(circularlagcorr(SPM.xX.X(4:(round(tn/150))+5,2),SPM.xX.X(4:(round(tn/150))+5,4),lags));
 plot(lags*1.5,r,'k','linewidth',2);
 [~,lagi]=max(r);
 peakLag=1.5*(lags(lagi));
line([peakLag peakLag],get(gca,'ylim'),'linestyle',':','color','k');
set(gca,'fontsize',14);
grid on
 

