function roi_tr_bined_pattern_granger_SL_perm(lagi)

disp(lagi)

loc='cluster';
set_parameters;
exp='merlin';
timeUnit='tr' ;
froidir='mor';
% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=0;
binSize_tested=[30]; % tr;
binStep=1;
lags_tested={-10:-1,-10:-4};
lags=lags_tested{lagi};

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');

for binSizei=1:length(binSize_tested);
    binSize=binSize_tested(binSizei);
    
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    for ri=1:length(rnames);
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0])-binSize+1);
        binN=length(1:binStep:(length(keptT)-binSize+1));
        voxn=size(data,1);
        
        if voxn*binSize>(length(lags)*2+1) ;
            for s=1:size(gdata,3);
                data_perm=gdata(:,:,s);
                gdata_perm=gdata;
                gdata_perm(:,:,s)=data;
                g=nanmean(gdata_perm,3);
                
                for bi=1:binN;
                    keptTi=((bi-1)*binStep+1):((bi-1)*binStep+binSize);
                    
                    y=g;
                    y=y(:,keptT(keptTi));
                    y=y(:);
                    
                    for li=1:length(lags);
                        Xl(:,:,li)=g(:,keptT(keptTi)+lags(li));
                        Xs(:,:,li)=data_perm(:,keptT(keptTi)+lags(li));
                    end
                    
                    Xl=reshape(Xl,roi_voxn*length(keptTi),length(lags));
                    Xs=reshape(Xs,roi_voxn*length(keptTi),length(lags));
                    
                    % centralized X
                    Xl=Xl-mean(Xl);
                    Xs=Xs-mean(Xs);
                    
                    [b_l_perm(ri,keptT(keptTi(1)),:,s),~,r_l,~,stats_l]=regress(y,[ones(size(Xl,1),1) Xl]);
                    [b_ls_perm(ri,keptT(keptTi(1)),:,s),~,r_ls,~,stats_ls]=regress(y,[ones(size(Xl,1),1) Xl Xs]);
                    [b_s_perm(ri,keptT(keptTi(1)),:,s),~,r_s,~,stats_s]=regress(y,[ones(size(Xl,1),1) Xs]);
                    
                    r2_l_perm(ri,keptT(keptTi(1)),s)=stats_l(1);
                    F_l_perm(ri,keptT(keptTi(1)),s)=stats_l(2);
                    p_l_perm(ri,keptT(keptTi(1)),s)=stats_l(3);
                    
                    r2_ls_perm(ri,keptT(keptTi(1)),s)=stats_ls(1);
                    F_ls_perm(ri,keptT(keptTi(1)),s)=stats_ls(2);
                    p_ls_perm(ri,keptT(keptTi(1)),s)=stats_ls(3);
                    
                    r2_s_perm(ri,keptT(keptTi(1)),s)=stats_s(1);
                    F_s_perm(ri,keptT(keptTi(1)),s)=stats_s(2);
                    p_s_perm(ri,keptT(keptTi(1)),s)=stats_s(3);
                    
                    
                    %The numerator of the F-statistic
                    rss_l = r_l'*r_l;
                    rss_s = r_s'*r_s;
                    rss_ls = r_ls'*r_ls;
                    
                    F_num = (rss_l- rss_ls)/length(lags);
                    %The denominator of the F-statistic
                    F_den = rss_ls/(length(y)-2*length(lags)-1);
                    %The F-Statistic
                    F_s2l_perm(ri,bi,s) = F_num/F_den;
                    
                    F_num = (rss_s- rss_ls)/length(lags);
                    F_l2l_perm(ri,bi,s) = F_num/F_den;
                    
                    clear Xl Xs
                end
            end
            
        else
            b_l_perm(ri,:,:,s)=NaN;
            b_ls_perm(ri,:,:,s)=NaN;
            b_s_perm(ri,:,:,s)=NaN;
            F_ls_perm(ri,:,s)=NaN;
            F_s_perm(ri,:,s)=NaN;
            F_l_perm(ri,:,s)=NaN;
            F_s2l_perm(ri,:,s)=NaN;
            F_l2l_perm(ri,:,s)=NaN;
            r2_ls_perm(ri,:,s)=NaN;
            r2_l_perm(ri,:,s)=NaN;
            r2_s_perm(ri,:,s)=NaN;
        end
    end
    save([expdir '/' exp '/fmri/pattern_regression_bined/' timeUnit '/roi/' froidir '/granger_SL_binsize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'b_l_perm','b_ls_perm','b_s_perm','F_l_perm','F_ls_perm','F_s2l_perm','F_l2l_perm','F_s_perm','r2_ls_perm','r2_l_perm','r2_s_perm','keptT','rnames');
    clear b F p r2 rnames b_l b_ls b_s r2_l r2_s r2_ls F_s2l F_l2l
end

