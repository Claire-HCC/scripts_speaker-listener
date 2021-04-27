clear all
loc='mypc'
set_parameters;
u=3; % 2:4; 3000 words are often not enough to create level 6 with more than to 2 consituents.
var=0.5;
speechR=1;
NRF={'linear'}; % {'linear', 'log','expn', 'triangle','time'};
pauseLen=3;
pauseEffect=0.1; % activity at pause period = min acitivity-sd of activity *pauseEffect
tb=readtable([expdir '/sherlock/sound/transcripts_sherlock_Claire.xlsx']);
dur_wd=round((tb.tmax-tb.tmin)*1000/speechR);
% speech R also applies to boundary
% duration
dur_boundary_sample=round(normrnd(pauseLen*1000,(pauseLen*1000)/6,1000,1));
dur_boundary_sample(dur_boundary_sample<0)=[];
% parameters of the associated normal distribution
mu = log((u^2)/sqrt(var+u^2));
sigma = sqrt(log(var/(u^2)+1));
constituentL=round(lognrnd(mu,sigma,1000,1));
constituentL(constituentL<=0)=[];

fsize=[19 20]
figure('unit','centimeter','position',[0 0 fsize])'

subplot(2,2,1);
histogram(dur_wd/1000,'facecolor',[0.5 0.5 0.5],'normalization','probability')
title({'Simulated word duration','based on Sherlock dataset'})
xlabel('word duration (sec)')
set(gca,'fontsize',14)
xlabel('Word duration (sec)')
ylabel('Probability (%)')
ylim([0 0.18])


subplot(2,2,3)
histogram(constituentL,'facecolor',[0.5 0.5 0.5],'normalization','probability');
ylabel('Probability (%)');
xlabel('Constituent length');
title({'Constituent length','mean = 3; variance = 0.5'})
set(gca,'fontsize',14)
set(gca,'xtick',1:7)
xlim([0.5 7.5])

subplot(2,2,4)
histogram(dur_boundary_sample/1000,'facecolor',[0.5 0.5 0.5],'normalization','probability')
ylabel('Probability (%)')
xlabel('Pause length (sec)')
title({'Pause length', 'mean = 3; standard deviation = 0.5'})
set(gca,'fontsize',14)


