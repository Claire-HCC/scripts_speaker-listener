clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
roi_ids=cell2mat(roi_table.id);
lags=-10:10;


for ei=3;%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/temporal_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bLagArea.mat'  ],'rnames','blag_pos','blag_neg');
    blag_temporal=blag_neg;

    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/SLeach/regression_SLeach_lag' num2str(min(lags)) '-' num2str(max(lags)) '_bLagArea.mat'  ],'rnames','blag_pos','blag_neg');
    blag_pattern=blag_neg;
    
    figure;

    for ri=1:length(rnames);
        text(blag_temporal(ri),blag_pattern(ri),rnames{ri},'fontsize',12);
    end
    xlabel({'Weighted average lag','temporal pattern'});
    ylabel({'Weighted average lag','temporal-spatial pattern'});
        xlim([-13 3]);
    ylim([-13 3]);
    set(gca,'xtick',-10:0);
    set(gca,'ytick',-10:0);
    title(exp)
end