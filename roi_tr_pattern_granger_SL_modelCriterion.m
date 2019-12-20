clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
crop_start=10;
lags=-[1:60];

for ei=[1 2 4];%1:2;
    exp=experiments{ei};
    
    rnames=table2array(roi_table(:,3));
    ris=find(cellfun(@(x) exist([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' x '.mat' ]),rnames)>0);
    rnames=rnames(ris);
    
    for ri=1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        keptT_s=find(([1:roi_tn]+min(lags))==1)+crop_start;
        keptT_e=min(roi_tn,roi_tn-max(lags));
        keptT=keptT_s:keptT_e;
        
        g=nanmean(gdata,3);
        y=g;
        y=y(:,keptT);
        y=y(:);
        
        for li=1:length(lags);
            % lagmatrix
            Xs(:,:,li)=data(:,keptT+lags(li));
            Xl(:,:,li)=g(:,keptT+lags(li));
            
            
            Xl_reshape=reshape(Xl,roi_voxn*length(keptT),li);
            Xs_reshape=reshape(Xs,roi_voxn*length(keptT),li);
            
            md_l=fitlm([ones(length(y),1) Xl_reshape],y);
            md_ls=fitlm([ones(length(y),1) Xl_reshape Xs_reshape],y);
            md_s=fitlm([ones(length(y),1) Xs_reshape],y);
            
            bic_ls(ri,li)= md_ls.ModelCriterion.BIC;
            bic_l(ri,li)= md_l.ModelCriterion.BIC;
            bic_s(ri,li)= md_s.ModelCriterion.BIC;
            
            aic_ls(ri,li)= md_ls.ModelCriterion.AIC;
            aic_l(ri,li)= md_l.ModelCriterion.AIC;
            aic_s(ri,li)= md_s.ModelCriterion.AIC;
            
            r_l=md_l.Residuals.Raw;
            r_s=md_s.Residuals.Raw;
            r_ls=md_ls.Residuals.Raw;
            
            rss_l = r_l'*r_l;
            rss_s = r_s'*r_s;
            rss_ls = r_ls'*r_ls;
            
            F_num = (rss_l- rss_ls)/length(lags);
            %The denominator of the F-statistic
            F_den = rss_ls/(length(y)-2*length(lags)-1);
            %The F-Statistic
            F_s2l(ri,li) = F_num/F_den;
            
         
            
        end
        clear Xl Xs
    end
    
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/granger_SL_modelCriterion_lag' num2str(min(lags)) '-' num2str(max(lags)) ''],'F_s2l','lags','rnames','keptT',...
        'aic_ls','aic_l','aic_s','bic_ls','bic_s','bic_l');
    clear F_s2l aic_ls aic_l aic_s bic_ls bic_s bic_l
end

toc
