clear all;
% F-test: https://support.sas.com/rnd/app/ets/examples/granger/index.htm
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
iters=1000;
lags=-10:-1; % how to determine lags?

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

for ei=3;%1:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    load([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'F_s2l','F_l2l');
    
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
    
    ris=find(ismember(rnames,'vPCUN'))
    for i=1:length(ris);%1:length(rnames);
        ri=ris(i);
        rname=rnames{ri};
        
        eventLabels_LG=h5read([expdir   exp '/fmri/hmm/' rname  '_findListenersEventInSpeaker.hdf5'],'///eventLabels_LG');
        K=length(unique(eventLabels_LG(:,1)));
        eventLabels_LG(1:(1-min(lags)))=0;
        eventLabels_LG_kept=unique(eventLabels_LG(eventLabels_LG~=0));
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        for iter=1:iters;
            if i==1;
                [data_perm, rp(:,iter)]=phase_rand2(data',1);
            else;
                data_perm=phase_rand2(data',1,rp(:,iter));
            end
            data_perm=data_perm';
            for s=1:size(gdata,3);
                gdata_perm(:,:,s)=phase_rand2(gdata(:,:,s)',1,rp(:,iter))';
            end
            
            g=nanmean(gdata_perm,3);
            for evi=1:K;
                keptT=find(eventLabels_LG==evi);
                
                if sum(keptT)>0;
                    y=g;
                    y=y(:,keptT);
                    y=y(:);
                    
                    for li=1:length(lags);
                        Xl(:,:,li)=g(:,keptT+lags(li));
                        Xs(:,:,li)=data_perm(:,keptT+lags(li));
                    end
                    
                    Xl=reshape(Xl,roi_voxn*length(keptT),length(lags));
                    Xs=reshape(Xs,roi_voxn*length(keptT),length(lags));
                    
                    % centralized X
                    Xl=Xl-mean(Xl);
                    Xs=Xs-mean(Xs);
                    
                    % add an constant
                    
                    [b_l_perm(ri,evi,:,iter),~,r_l,~,stats_l]=regress(y,[ones(size(Xl,1),1)  Xl]);
                    [b_ls_perm(ri,evi,:),~,r_ls,~,stats_ls]=regress(y,[ones(size(Xl,1),1) Xl Xs]);
                    [b_s_perm(ri,evi,:,iter),~,r_s,~,stats_s]=regress(y,[ones(size(Xl,1),1) Xs]);
                    
                    r2_l_perm(ri,evi,iter)=stats_l(1);
                    F_l_perm(ri,evi,iter)=stats_l(2);
                    p_l_perm(ri,evi,iter)=stats_l(3);
                    
                    r2_ls_perm(ri,evi,iter)=stats_ls(1);
                    F_ls_perm(ri,evi,iter)=stats_ls(2);
                    p_ls_perm(ri,evi)=stats_ls(3);
                    
                    r2_s_perm(ri,evi,iter)=stats_s(1);
                    F_s_perm(ri,evi,iter)=stats_s(2);
                    p_s_perm(ri,evi,iter)=stats_s(3);
                    
                    %The numerator of the F-statistic
                    rss_l = r_l'*r_l;
                    rss_s = r_s'*r_s;
                    rss_ls = r_ls'*r_ls;
                    
                    F_num = (rss_l- rss_ls)/length(lags);
                    %The denominator of the F-statistic
                    F_den = rss_ls/(length(y)-2*length(lags)-1);
                    %The F-Statistic
                    F_s2l_perm(ri,evi,iter) = F_num/F_den;
                    
                    F_num = (rss_s- rss_ls)/length(lags);
                    %The denominator of the F-statistic
                    F_den = rss_ls/(length(y)-2*length(lags)-1);
                    %The F-Statistic
                    F_l2l_perm(ri,evi,iter) = F_num/F_den;
                    
                    clear Xl Xs
                end
            end
            
        end
        p_s2l(ri,:)=sum(squeeze(F_s2l_perm(ri,:,:))'>F_s2l(ri,:))/iters;
        p_l2l(ri,:)=sum(squeeze(F_l2l_perm(ri,:,:))'>F_l2l(ri,:))/iters;
    end
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_event_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm'],'b_l_perm','b_ls_perm','b_s_perm','F_l_perm','r2_l_perm','p_ls_perm','F_ls_perm','r2_ls_perm','p_ls_perm','F_s_perm','r2_s_perm','p_s_perm','F_s2l_perm','eventLabels_LG_kept','rnames','p_s2l','p_l2l');
    % clear b F p r2 rnames coupling
    
end
toc
beep
