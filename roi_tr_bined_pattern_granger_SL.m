function roi_tr_bined_pattern_granger_SL(lagi)

disp(lagi)

loc='cluster';
set_parameters;
exp='merlin';
timeUnit='tr' ;
froidir='mor';
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
binSize_tested=[ 15 30]; % tr;
lags_tested={-10:-1,-10:-4};
lags=lags_tested{lagi};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
rnames=table2array(roi_table(:,3));
ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
rnames=rnames(ris);

for binSizei=1:length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
            load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
            g=nanmean(gdata,3);
            
            roi_voxn=size(gdata,1);
            tn=size(gdata,2);
            
            for t=1:tn;
                t_bin=t:(t+binSize-1);
                
                if min(t_bin)+min(lags)>=1 & t_bin+max(lags)<=tn & max(t_bin)<=tn;
                    y=g;
                    y=y(:,t_bin);
                    y=y(:);
                    
                    for li=1:length(lags);
                        Xl(:,:,li)=g(:,t_bin+lags(li));
                        Xs(:,:,li)=data(:,t_bin+lags(li));
                    end
                    
                    Xl=reshape(Xl,roi_voxn*length(t_bin),length(lags));
                    Xs=reshape(Xs,roi_voxn*length(t_bin),length(lags));
                    
                    % centralized X
                    Xl=Xl-mean(Xl);
                    Xs=Xs-mean(Xs);
                    
                    [b_l(ri,t,1:(1+length(lags))),~,r_l,~,stats_l]=regress(y,[ones(size(Xl,1),1) Xl]);
                    [b_ls(ri,t,1:(1+length(lags)*2)),~,r_ls,~,stats_ls]=regress(y,[ones(size(Xl,1),1) Xl Xs]);
                    [b_s(ri,t,1:(1+length(lags))),~,r_s,~,stats_s]=regress(y,[ones(size(Xl,1),1) Xs]);
                    
                    r2_l(ri,t)=stats_l(1);
                    F_l(ri,t)=stats_l(2);
                    p_l(ri,t)=stats_l(3);
                    
                    r2_ls(ri,t)=stats_ls(1);
                    F_ls(ri,t)=stats_ls(2);
                    p_ls(ri,t)=stats_ls(3);
                    
                    r2_s(ri,t)=stats_s(1);
                    F_s(ri,t)=stats_s(2);
                    p_s(ri,t)=stats_s(3);
                    
                    
                    %The numerator of the F-statistic
                    rss_l = r_l'*r_l;
                    rss_s = r_s'*r_s;
                    rss_ls = r_ls'*r_ls;
                    
                    F_num = (rss_l- rss_ls)/length(lags);
                    %The denominator of the F-statistic
                    F_den = rss_ls/(length(y)-2*length(lags)-1);
                    %The F-Statistic
                    F_s2l(ri,t) = F_num/F_den;
                    
                    F_num = (rss_s- rss_ls)/length(lags);
                    F_l2l(ri,t) = F_num/F_den;
                    
                    clear Xl Xs
                    
                else
                    b_l(ri,t,1:(1+length(lags)))=NaN;
                    b_ls(ri,t,1:(1+length(lags)*2))=NaN;
                    b_s(ri,t,1:(1+length(lags)))=NaN;
                    F_ls(ri,t)=NaN;
                    F_s(ri,t)=NaN;
                    F_l(ri,t)=NaN;
                    F_s2l(ri,t)=NaN;
                    F_l2l(ri,t)=NaN;
                    r2_ls(ri,t)=NaN;
                    r2_l(ri,t)=NaN;
                    r2_s(ri,t)=NaN;
                end
            end
            disp(ri)
    end
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'b_l','b_ls','b_s','F_l','F_ls','F_s2l','F_l2l','F_s','r2_ls','r2_l','r2_s','rnames');
    clear b F p r2 b_l b_ls b_s r2_l r2_s r2_ls F_s2l F_l2l
end

