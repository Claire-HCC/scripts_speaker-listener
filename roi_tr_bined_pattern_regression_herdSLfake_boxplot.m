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

rnames_selected={'HG_L','vPCUN','aANG_R','aCUN'};
figure;

for rii=1:length(rnames_selected);
    rname_selected=rnames_selected{rii};
    herdz_ds=nan(4,48);
    for ei=1:4;%1:4;
        exp=experiments{ei};
        
        load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'herdz_d','rnames');
        ri=find(ismember(rnames,rname_selected));
        rname=rnames{ri};
        
        
        herdz_ds(ei,1:size(herdz_d,2))=herdz_d(ri,:);
        
    end
    subplot(2,2,rii)
    boxplot(herdz_ds');
    ylim([-1 1]);
  
    set(gca,'xticklabels',experiments);
    ylabel('Real > Pseudo Herding efect(?)');
    title([rname]);
    hold on
    line([0 5],[0 0 ],'Color','k');
end

