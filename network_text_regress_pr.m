clear all
loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='restFc_isc30PercMasked_75Overlap_cluster6';
networks={'AUD','vLAN','dLAN','DMNa','Attention','DMNb'};
win=-5:15;

crop_start=25;
crop_end=20;

% xBF.dt=1.5;
% xBF.name='hrf (with time and dispersion derivatives)';
% bf = spm_get_bf(xBF);
% bf=bf.bf;
bf=eye(length(win));

for ei=6;
    exp=exp_parameters.experiments{ei};
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' networks{1} '.mat' ],'gdata','keptvox');
    % keep only voxel with peakLag 0
    [~,tn,listenerN]=size(gdata);
    keptT=(crop_start+1):(tn-crop_end);
    
    load([expdir exp '/sound/onsets.mat'],'offsets_tr_pr');
    
    load([expdir exp '/sound/' exp '_listener_audenv.mat'],'aud')
    aud=zscore(aud);
    
    X_temp=zeros(tn,3);
    X_temp(round(offsets_tr_pr),1)=1;
    
    X_temp=X_temp((abs(min(win(win<0)))+1):end,:);
    X_temp=X_temp-mean(X_temp);
    
    for xi=1:size(X_temp,2);
        for bi=1:size(bf,2);
            if xi==1 & bi==1;
                X=conv(X_temp(:,xi),bf(:,bi));
            else
                X(:,(xi-1)*size(bf,2)+bi)=conv(X_temp(:,xi),bf(:,bi));
            end
        end
    end
    
    X=[aud(1:tn)  X(1:tn,:)];
    X=X-nanmean(X);
    X=[ones(tn,1) X];
    
    for sdi=1:length(networks);
        seed=networks{sdi};
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '/listenerAll_' seed '.mat' ],'gdata','keptvox');
        % keep only voxel with peakLag 0
        
        [~,tn,listenerN]=size(gdata);
        keptT=(crop_start+1):(tn-crop_end);
        
        gdata(:,~ismember(1:tn,keptT),:)=NaN;
        gdata(:,:,exp_parameters.subjects_excluded{ei})=NaN;
        gdata=nanmean(gdata,1);
        resid=[];
        b=[];
        Y=[];
        for si=1:listenerN;
            y=squeeze(gdata(:,:,si))';
            [b,bint,~,~,stats] = regress(y,X);
            Y(1,:,si)=X(:,2:end)*b(2:end);
            % keep the intercept
            resid(1,:,si)=y-X(:,2:end)*b(2:end);
        end
        y=gdata;
        gdata=resid;
        
        mkdir([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '_audPrpauseResid/']);
        save([expdir '/' exp '/fmri/timeseries/' timeUnit '/network/' froidir '_audPrpauseResid/listenerAll_' seed '.mat' ],'gdata','keptvox','Y','y');
    end
end




