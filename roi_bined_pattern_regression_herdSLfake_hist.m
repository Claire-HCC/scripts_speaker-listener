% bined r2 vs std between subjects
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

% -4 for merlin, -3 for sherlock
lags=-7:-4;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
fsize=[34 30];
fig1=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
fig2=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
fig3=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);


rnames_selected={'HG_L','vPCUN','aCUN'};


% [~,ri]=max(herdm);

for ei=1:4;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herdz_d','lags','herd','herd_p','herd_sig_fdr_pos','herd_sig_fdr_neg','herd_null','rnames');
    ris=find(ismember(rnames,rnames_selected));
    
    for rii=1:length(ris);
        ri=ris(rii);
        rname=rnames{ri};
        eval(['figure(fig' num2str(rii) ')']);
        
        herdm=mean(herd,2);
        herdm_null=nanmean(herd_null,2);
        rnames(herd_sig_fdr_pos==1)
        subplot(2,2,ei);
        %         histogram(herd(ri,:),'BinWidth',0.05)
        %         hold on
        %         histogram(herd_null(ri,:),'BinWidth',0.05);
        histogram(herdz_d(ri,:),'BinWidth',0.05)
        xlim([-0.8 0.8]);
        hold off
        title([exp ' ' rname]);
        grid on
    end
end
legend('real herding effect','pseud-herding effect');