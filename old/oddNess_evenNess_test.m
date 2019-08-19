clear all
% close all
loc='mypc';
set_parameters

load(['Y:\claire\speaker-listener\merlin\sound\listenerEvents.mat']);

rname='vPCUN';
role='listener';
exp='merlin';
load([expdir exp '/fmri/mat/roi/segment/segment_' role '_' rname '.mat']);

data=zscore(data,0,2);
for s=1:18;
    train_subjects=1:18;
    train_subjects=train_subjects(train_subjects~=s);
    train_data=data(:,:,train_subjects);
    test_data=data(:,:,s);
    
    for segi=1:22;
        ai=1:2:21;
        bi=2:2:22;
        
        train_A(:,segi,s)=mean(mean(train_data(:,ai(ai~=segi),:),2),3);
        train_B(:,segi,s)=mean(mean(train_data(:,bi(bi~=segi),:),2),3);
    end
    
    test_A(1,:,s)=corr_col(test_data,train_A(:,:,s));
    test_B(1,:,s)=corr_col(test_data,train_B(:,:,s));
end

for iter=1:100;
    
    [~,~,segment_shuffled]=unique(rand(22,1));
    for s=1:18;
        train_subjects=1:18;
        
        train_subjects=train_subjects(train_subjects~=s);
        test_data=data(:,:,s);
        
        for segi=1:22;
            train_data=data(:,:,train_subjects);
            train_data(:,segi,:)=NaN;
            train_data=train_data(:,segment_shuffled,:);
            
            ai=1:2:21;
            bi=2:2:22;
            
            train_A(:,segi,s)=nanmean(nanmean(train_data(:,ai(ai~=segi),:),2),3);
            train_B(:,segi,s)=nanmean(nanmean(train_data(:,bi(bi~=segi),:),2),3);
        end
        
        test_A_perm(1,:,s,iter)=corr_col(test_data,train_A(:,:,s));
        test_B_perm(1,:,s,iter)=corr_col(test_data,train_B(:,:,s));
    end
end

save(['oddNess_evenNess_' rname '.mat'],'test_A_perm','test_A_perm');

test_A_perm_m=squeeze(mean(test_A_perm,3));
test_B_perm_m=squeeze(mean(test_B_perm,3));

test_A_m=squeeze(mean(test_A,3));
test_B_m=squeeze(mean(test_B,3));

test_A_ci95=quantile(test_A_perm_m',[ (1-0.95)/2 1-(1-0.95)/2]);
test_B_ci95=quantile(test_B_perm_m',[ (1-0.95)/2 1-(1-0.95)/2]);


figure;
subplot(1,2,1)
hold on
ciplot(test_A_ci95(1,ai),test_A_ci95(2,ai),1:length(ai),[1 0.1 0.1],0.1);
p1=plot(test_A_m(ai));
ciplot(test_B_ci95(1,ai),test_B_ci95(2,ai),1:length(ai),[0.1 0.1 1],0.1);
p1=plot(test_B_m(ai));

hold off
grid on
ylim([-1 1]);
xlim([1 11]);
set(gca,'xtick',1:11);
set(gca,'xticklabel',ai);
xlabel('test odd');

subplot(1,2,2)
ciplot(test_A_ci95(1,bi),test_A_ci95(2,bi),1:length(bi),[1 0.1 0.1],0.1);
hold on
p1=plot(test_A_m(bi));
ciplot(test_B_ci95(1,bi),test_B_ci95(2,bi),1:length(bi),[0.1 0.1 1],0.1);
p2=plot(test_B_m(bi))
hold off
grid on
ylim([-1 1]);
xlim([1 11]);
xlabel('test even');
legend([p1 p2],'train odd','train even')
set(gca,'xtick',1:11);
set(gca,'xticklabel',bi);



figure;
subplot(1,2,1)
hold on

ciplot(test_A_ci95(1,:),test_A_ci95(2,:),1:size(test_A_ci95,2),[1 0.1 0.1],0.1);
p1=plot(test_A_m);
ciplot(test_B_ci95(1,:),test_B_ci95(2,:),1:size(test_A_ci95,2),[0.1 0.1 1],0.1);
p2=plot(test_B_m);

hold off
grid on
ylim([-1 1]);
set(gca,'xtick',1:22);
set(gca,'xticklabel',1:22);
xlabel('test');
legend([p1 p2],'train odd','train even')


subplot(1,2,2)
d_corr=corr(zscore(mean(data,3),0,2));
d_corr(eye(22)==1)=NaN;
boxplot(d_corr,'PlotStyle','compact','color','k')
grid on;
% ylim([-0.5 0.5]);
xlim([0 22]);
xlabel('segment')
ylabel('correlation with other segments')
set(gca,'xtick',1:22);
set(gca,'xticklabel',1:22);


% acc graph
labels=repmat([1 2],1, 11,18);
ab=((test_A-test_B));
labels_predicted=zeros(size(labels));
labels_predicted(ab>0)=1;
labels_predicted(ab<0)=2;
acc=(labels_predicted(:,1:22,:)==labels(:,1:22,:));
% m=100*squeeze(mean(acc,3));
% ci95=100*ci(squeeze(acc)',0.95);
[m,ci95] = binofit(squeeze(sum(acc,3)),18,0.05)
m=m*100;
ci95=100*[ci95(:,2)-ci95(:,1)];

figure
subplot(1,2,1);
bar(m,'k');
hold on
errorbar(m,ci95,'.','color','r');
hold off
line([0  30],[50 50],'color',[0.7 0.7 0.7])
xlim([0 31]);
xlabel('segment');
ylim([0 100]);
ylabel('Acc (%)');
title({rname,'odd/even part', ['mean ACC: ' num2str(round(mean(m))) '%']});


subplot(1,2,2);
labels=repmat([1 2],1, 11,18,iter);
ab=((test_A_perm-test_B_perm));
labels_predicted=zeros(size(labels));
labels_predicted(ab>0)=1;
labels_predicted(ab<0)=2;
acc=(labels_predicted==labels);
m=100*(squeeze(mean(mean(acc,3),4)));
ci95=100*ci(reshape(squeeze(acc),22,18*100)',0.95);
bar(m,'k');
hold on
errorbar(m,ci95,'.','color','r');
hold off
line([0  30],[50 50],'color',[0.7 0.7 0.7])
xlim([0 31]);
xlabel('segment');
ylim([0 100]);
ylabel('Acc (%)');
title({rname,'odd/even part', ['mean ACC: ' num2str(round(mean(m))) '%']});




