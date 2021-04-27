%% ttest across subject
clear all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags_tested={-10:10,  -40:40};

rname='pANG_L';
fsize=[30,8];
figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
for ei=1:4;
    exp=experiments{ei};
    
    for lagi=1%:length(lags_tested);
        lags=lags_tested{lagi};
        
        ri=find(ismember(rnames,rname));
        
        load([expdir '/' exp '/fmri/temporal/lagcorr/' timeUnit '/roi2roi/' froidir '/LL_leave1out/lag' num2str(min(lags)) '-' num2str(max(lags)) ],'rnames','r','lags','keptT');

        subplot(1,4,ei);
        plot(lags,squeeze(r(ri,:,:))','k');
        xlabel('lags (TR)');
        ylabel('R');
        title({rname,exp});
        set(gca,'fontsize',14);
ylim([-0.4 0.8]);
grid on
    end
end