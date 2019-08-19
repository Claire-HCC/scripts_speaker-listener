clear all;
close all
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';

close all; for ci=1:7;
figure; plot(lags,mean(b(clusters==ci,2:end),1));
end

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
lags=-40:40;
smoothSpan=3;
for ei=1:2;
    exp=experiments{ei};
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F','b','rnames');
    load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag'  num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'F_perm','b_perm','rnames');
    
    
    p=sum(F<F_perm,3)/1000;
    sig=(p<(0.05/length(rnames)));
    sig=fdr0(p,0.05);
    
    networks=unique(roi_table.network);
    for ni=1:length(networks);
        %         if ei==1; f(ni)=figure; else ;set(0, 'currentfigure', f(ni));  end
        %         hold on
        network=networks{ni};
        b_temp(ni,:,ei)=nanmean(b(sig==1 & ismember(rnames, roi_table.region(ismember(roi_table.network,network))),2:end ));
        
        if ei==2 & length(b_temp)>0;
            figure
            % b_temp=b( ismember(rnames,network),2:end );
            r=corr(squeeze(b_temp(ni,:,:)));
            b_m=zscore(mean(b_temp(ni,:,:),3));
            y1=smooth(zscore(squeeze(b_temp(ni,:,1)),0,2),smoothSpan);
            plot(lags*tr(ei),y1);
            hold on 
            y2=smooth(zscore(squeeze(b_temp(ni,:,2)),0,2),smoothSpan);
             plot(lags*tr(ei),y2);
            hold off
            title([network '; r =' sprintf('%02f',r(2,1))]);
            grid on
            xlabel('speaker precedes    shift(sec)    listern precedes');
            xlim([-40 40])
            % hold off
            legend(experiments(1:2));
        end
    end
end
