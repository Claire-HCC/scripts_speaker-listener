clear all;
% F-test: https://support.sas.com/rnd/app/ets/examples/granger/index.htm
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data

lags_tested={-10:-1,-10:-4,-60:-4};


load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

for li=1:length(lags_tested);
    lags=lags_tested{li};
    
    for ei=3;%1:4;
        exp=experiments{ei};
        rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
        rnames={rnames.name};
        rnames=strrep(rnames,'zscore_listenerAll_','');
        rnames=strrep(rnames,'.mat','');
        rnames=rnames';
        rnames=rnames(ismember(rnames,roi_table.region));
        
        roi_ids=[];
        
        b_l=[];
        r2_l=[];
        F_l=[];
        b_ls=[];
        r2_ls=[];
        F_ls=[];
        
        b_s=[];
        r2_s=[];
        F_s=[];
        
        % ris=find(ismember(rnames,'vPCUN'))
        % for i=1:length(ris);%1:length(rnames);
        %  ri=ris(i);
        
        for ri=1:length(rnames);
            rname=rnames{ri};
            
            eventLabels_LG=h5read([expdir   exp '/fmri/hmm/findListenersEventInSpeaker/' rname  '.hdf5'],'///eventLabels_LG');
            K=length(unique(eventLabels_LG(:,1)));
            eventLabels_LG(1:(1-min(lags)))=0;
            eventLabels_LG_kept{ri}=unique(eventLabels_LG(eventLabels_LG~=0));
            
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
            
            roi_voxn=size(gdata,1);
            roi_tn=size(gdata,2);
            
            g=nanmean(gdata,3);
            
            b_l{ri}=NaN(max(eventLabels_LG_kept{ri}),length(lags)+1);
            b_ls{ri}=NaN(max(eventLabels_LG_kept{ri}),length(lags)*2+1);
            b_s{ri}=NaN(max(eventLabels_LG_kept{ri}),length(lags)+1);
            F_ls{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            F_s{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            F_l{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            F_s2l{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            F_l2l{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            r2_ls{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            r2_l{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            r2_s{ri}=NaN(max(eventLabels_LG_kept{ri}),1);
            
            for evi=1:K;
                keptT=find(eventLabels_LG==evi);
                
                if voxn*length(keptT)>(length(lags)*2+1) ;
                    
                    if sum(keptT)>0;
                        y=g;
                        y=y(:,keptT);
                        y=y(:);
                        
                        for li=1:length(lags);
                            Xl(:,:,li)=g(:,keptT+lags(li));
                            Xs(:,:,li)=data(:,keptT+lags(li));
                        end
                        
                        Xl=reshape(Xl,roi_voxn*length(keptT),length(lags));
                        Xs=reshape(Xs,roi_voxn*length(keptT),length(lags));
                        
                        % centralized X
                        Xl=Xl-mean(Xl);
                        Xs=Xs-mean(Xs);
                        
                        % add an constant
                        
                        [b_l{ri}(evi,:),~,r_l,~,stats_l]=regress(y,[ones(size(Xl,1),1)  Xl]);
                        [b_ls{ri}(evi,:),~,r_ls,~,stats_ls]=regress(y,[ones(size(Xl,1),1) Xl Xs]);
                        [b_s{ri}(evi,:),~,r_s,~,stats_s]=regress(y,[ones(size(Xl,1),1) Xs]);
                        
                        r2_l{ri}(evi)=stats_l(1);
                        F_l{ri}(evi)=stats_l(2);
                        p_l{ri}(evi)=stats_l(3);
                        
                        r2_ls{ri}(evi)=stats_ls(1);
                        F_ls{ri}(evi)=stats_ls(2);
                        p_ls{ri}(evi)=stats_ls(3);
                        
                        r2_s{ri}(evi)=stats_s(1);
                        F_s{ri}(evi)=stats_s(2);
                        p_s{ri}(evi)=stats_s(3);
                        
                        %The numerator of the F-statistic
                        rss_l = r_l'*r_l;
                        rss_s = r_s'*r_s;
                        rss_ls = r_ls'*r_ls;
                        
                        F_num = (rss_l- rss_ls)/length(lags);
                        %The denominator of the F-statistic
                        F_den = rss_ls/(length(y)-2*length(lags)-1);
                        %The F-Statistic
                        F_s2l{ri}(evi) = F_num/F_den;
                        
                        F_num = (rss_s- rss_ls)/length(lags);
                        F_l2l{ri}(evi) = F_num/F_den;
                        
                    end
                end
                clear Xl Xs
                
            end
        end
        save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_l','b_ls','b_s','F_l','F_ls','F_s2l','F_l2l','F_s','r2_ls','r2_l','r2_s','rnames','eventLabels_LG_kept');
        % clear b F p r2 rnames coupling
    end
end

toc
beep
