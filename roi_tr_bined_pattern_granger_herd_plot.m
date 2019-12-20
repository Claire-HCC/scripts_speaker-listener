% bined r2 vs std between subjects
% clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
rname='vPCUN'
% -4 for merlin, -3 for sherlock
lags_tested={-10:-1,-10:-4};
binSize_tested=[  10 15 30 ];
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');

% load([expdir '/scripts_speaker-listener/merlin_' rname '_hmm_findListenersEventInSpeaker.mat'], 'segments_LG');

%r2_l leads to lower r_ls-r2_l; which might not be fair to the speaker,
%because listner can predict, but not as good as speaker.
% the solutions are:
% 1. set lags=-10:-4
% 2. use F-value as coupling index
% don't use granger causality.

for ei=3%:4;%1:4;
    exp=experiments{ei};
    fsize=[20 9];
    % f1=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    % f2=figure('unit','centimeter','position',[0 0 fsize],'paperposition',[0 0 fsize],'papersize',fsize);
    
    for bi=1%:length(binSize_tested);
        binSize=binSize_tested(bi);
        binStep=1;%binSize;
        
        for lagi=2;%1:length(lags_tested);
            lags=lags_tested{lagi};
            
             save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/herd_granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats'],'p','herd_null','herd','rnames','lags','keptT','sig_fdr');
             p_herd=p;
             
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F_s2l','r2_l','r2_ls','rnames');
            coupling_sl=r2_ls-r2_l;
            
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats'],'p');
            p_sl=p;
            
            load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/regression_LL2_binSize' num2str(binSize) '_lag0-0' ],'F','r2');
            coupling_ll=nanmean(r2,3);
            
            keptT=find(sum(isnan(coupling_sl),1)==0  & sum(isnan(coupling_ll),1)==0 );
            keptT=keptT((binStep/2):binStep:end);
            
            load([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/granger_SL_lag' num2str(min(lags)) '-' num2str(max(lags)) '_stats'],'p');
            
            % ris=find(ismember(rnames,{'HG_L','vPCUN','aCUN'}));
            ris=find(ismember(rnames,rname));
            % [~,ri]=max(herdm);
            
            for rii=1:length(ris);
                ri=ris(rii);
                
                eventLabels_LG=h5read([expdir   exp '/fmri/hmm/findListenersEventInSpeaker/' rname  '.hdf5'],'///eventLabels_LG');
                K=length(unique(eventLabels_LG(:,1)));
                % cols=distinguishable_colors(K);
                
                coupling_sl_temp=coupling_sl(ri,keptT)';
                coupling_ll_temp=nanmean(coupling_ll(ri,keptT,:),3)';
                sig_sl_temp=p_sl(ri,keptT)==0;
                
                figure
                % set(0, 'currentfigure', f1);
                % subplot(length(binSize_tested),length(lags_tested),(bi-1)*3+lagi);
                hold on
                
                legend({'SL coupling','LL coupling'},'fontsize',13); legend boxoff
                
                cols=zeros(length(keptT),3);
                colx=sig_sl_temp;
                coly=coupling_ll_temp>((max(coupling_ll_temp)-min(coupling_ll_temp))/2+min(coupling_ll_temp));
                
                
                for ti=1:length(keptT);
                    %    coli=1+ceil(atan2d(coupling_sl_temp(ti)>median(coupling_sl_temp),coupling_ll_temp(ti)-median(coupling_ll_temp)));
                    if  colx(ti)==1 & coly(ti)==1;
                        cols(ti,:)=[1 0.7 0];
                    elseif colx(ti)==1 & coly(ti)==0;
                        cols(ti,:)=[0 1 0];
                    elseif colx(ti)==0 & coly(ti)==0;
                        cols(ti,:)=[0 0 1];
                    elseif colx(ti)==0 & coly(ti)>0;
                        cols(ti,:)=[0.7 0.1 1];
                    end
                    area((keptT(ti)-binStep/2+1):(keptT(ti)+binStep/2),ones(1,binStep),'Facecolor',cols(ti,:),'FaceAlpha',0.3,'linestyle','none');
                end
                p1=plot(keptT,[coupling_sl_temp],'r','linewidth',2);
                plot(keptT,[coupling_ll_temp],'k','linewidth',2);
                
                title({[exp ', ' rnames{ri}],['binSize' num2str(binSize) ', Lags:' num2str(min(lags)) '~' num2str(max(lags)) ]});
                xlim([min(keptT)-binSize/2 max(keptT)+binSize/2]);
                ylim([0 max(coupling_sl_temp(:))+0.05]);
                
                % ylim([0 1.1*max(max(squeeze(coupling_slf(ri,:,:))))]);
                xlabel('Time (TR)'); % ylabel('Coupling');
                set(gca,'fontsize',13);
                
                
                figure
                hold on
                %  set(0, 'currentfigure', f2);
                %        subplot(length(binSize_tested),length(lags_tested),(bi-1)*3+lagi);
                for ti=1:length(keptT);
                    coli=1+ceil(atan2d(coupling_sl_temp(ti)-median(coupling_sl_temp),coupling_ll_temp(ti)-median(coupling_ll_temp)));
                    if coli<0; coli=360+coli;end
                    scatter(coupling_sl_temp(ti),coupling_ll_temp(ti),60,cols(ti,:),'filled');
                end
                
                % scatter(coupling_sl_temp,coupling_ll_temp,60,'k','filled');hold on
                % scatter(coupling_sl_temp(p_temp==0),coupling_ll_temp(ri,p_temp==0)',60,'b','filled');
                [r ]=corr(coupling_sl_temp,coupling_ll_temp);
                xlabel('SL coupling'); ylabel('LL coupling');
                hold off
                % if herd_sig_fdr_pos(ri)==1; star='*'; else  star='';end
                % title(sprintf('Averaged R = %.02f%s',herdm(ri),star));
                
                title({sprintf('Averaged R = %.2f, p = %.3f',r,p_herd(ri)),['binSize' num2str(binSize) ', Lags:' num2str(min(lags)) '~' num2str(max(lags)) ]});
                set(gca,'fontsize',13);
            end
        end
    end
    clear herd herd_null herd_p herd_sig_fdr_pos herd_sig_fdr_neg
end

%y=cos(0:((2*pi)/359):2*pi);
%x=sin(0:((2*pi)/359):2*pi);
%scatter(x,y,40,cols)