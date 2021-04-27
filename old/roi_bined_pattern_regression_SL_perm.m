clear all;

tic
 loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;
binSize=30; % tr;
iters=1000;
lags=-4:4;
for ei=3%1:4;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/zscore_listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'zscore_listenerAll_','');
    rnames=strrep(rnames,'.mat','');
    rnames=rnames';
    rnames=rnames(ismember(rnames,roi_table.region));
    
    b=[];
    roi_ids=[];
    r2=[];
    F=[];
    for ri=1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_speaker_' rname '.mat'],'data');
        
        
        %% isc or pattern
        % gdata=nanmean(gdata,1);
        % data=nanmean(data,1);
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        
        g=nanmean(gdata,3);
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags, 0]));
        binN=floor(length(keptT)/binSize);
        keptT=keptT(1:(binSize*binN));
        for iter=1:iters;
            clear data_mat
            data_perm=zeros(size(data))';
            if ri==1;
                 [data_perm((crop_start+1):end,:), rp(:,iter)]=phase_rand2(data(:,(crop_start+1):end)',1);
            else;
                [data_perm((crop_start+1):end,:)]=phase_rand2(data(:,(crop_start+1):end)',1,rp(:,iter));
            end
            data_perm=data_perm';
            
            for bi=1:binN;
                keptTi=((bi-1)*binSize+1):(bi*binSize);
                
                % substract the mean pattern
                y=zscore(g(:,keptT(keptTi)),0,2);
                y=y(:);
                for li=1:length(lags);
                     % substract the mean pattern
                    X(:,:,li)=zscore(data_perm(:,keptT(keptTi)+lags(li)),0,2);
                end
                
                X=reshape(X,roi_voxn*length(keptTi),length(lags));
                
                % centralized X
                X=X-mean(X);
                
                % add an constant
                X=[ones(size(X,1),1) X];
                
                [b_perm(ri,bi,:,iter),~,~,~,stats]=regress(y,X);
                
                r2_perm(ri,bi,iter)=stats(1);
                F_perm(ri,bi,iter)=stats(2);
                p_perm(ri,bi,iter)=stats(3);
                
                Y=X*squeeze(b_perm(ri,bi,:,iter));
                Y=reshape(Y,roi_voxn,length(keptTi));
                y=reshape(y,roi_voxn,length(keptTi));
                %  coupling(ri,:)=zeros(1,length(keptT));
                coupling(ri,keptT(keptTi),iter)=corr_col(y,Y);
                
                clear X
            end
        end
    end
    couplingz_perm=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/bined_pattern_regression/' timeUnit '/roi/' froidir '/regression_SL_bin' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '_perm' ],'b_perm','F_perm','r2_perm','p_perm','couplingz_perm','lags','rnames','keptT','binSize','binN');
    clear b F p r2 rnames coupling r
end

toc
beep
