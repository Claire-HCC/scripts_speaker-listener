
clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-40:40;

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};
% load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/herd_exps_' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],'b_sl','b_ll','keptT','rnames','lags','herd','herd_p','herd_sig');
timeBinN=2;

ris=[3 54];
% ris=find(sum(herd_sig==0,2)==0);
for i=1:length(ris);
    ri=ris(i);
    
    fsize=[35 20];
    fig1=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    for ei=1:4;
        exp=experiments{ei};
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags))],'couplingz','keptT','rnames');
        cp_ll=nanmean(couplingz,3);
        cp_ll=cp_ll(:,keptT);
        
        load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags))],'couplingz','keptT');
        cp_sl=couplingz;
        cp_sl=cp_sl(:,keptT);
        
        timeBinLen=floor(size(cp_sl,2)/timeBinN);
        
        rname=rnames{ri};
        for timeBinInd=1:timeBinN;
            timei=(timeBinLen*(timeBinInd-1)+1):(timeBinLen*timeBinInd);
            
            figure(fig1);
            subplot(timeBinN,4,ei+(timeBinInd-1)*4);
            [r p ]=corr(cp_sl(ri,timei)',cp_ll(ri,timei)');
            
            p = polyfit(cp_sl(ri,timei)',cp_ll(ri,timei)', 1);
            Y = polyval(p, (cp_sl(ri,timei)'));
            
            scatter(cp_sl(ri,timei),cp_ll(ri,timei),20,'k','filled');
            hold on
            plot(cp_sl(ri,timei),Y,'color',[0.5 0.5 0.5],'linewidth',1.5);
            title([rname ', ' exp ', ' sprintf('R=%.02f',r)]);
            hold off
            if ei ==1;
                ylabel({['timeBin' num2str(timeBinInd)],'Speaker-Listerner Coupling'});
            else
                ylabel('Listener-Listerner Coupling');
            end
            xlabel('Speaker-Listerner Coupling');
            
        end
    end
end


