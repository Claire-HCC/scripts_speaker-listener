
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

crop_start=10;
lags=-150:150;
for ei=1:2;
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
    for ri=[3 54 57]%1:length(rnames);
        clear data_mat
        rname=rnames{ri};
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/zscore_listenerAll_' rname '.mat' ],'gdata');
        gdata=mean(gdata,1);
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));

        
        for s=1:size(gdata,3);
            
            othersi=1:size(gdata,3);
            others=othersi(othersi~=s);
            others=nanmean(gdata(:,:,othersi),3);
          %  others(:,(crop_start+1):end)=zscore( others(:,(crop_start+1):end),0,2);
            
            self=gdata(:,:,s);
          %  self(:,keptT)=zscore(self(:,keptT),0,2);
            
            y=others;
            y=y(:,keptT);
            y=y(:);
           
            
            for li=1:length(lags);
                X(:,:,li)=self(:,keptT+lags(li));
            end
           
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            X=X-mean(X);
            X=[ones(size(X,1),1) X];
            [b(ri,:,s),~,~,~,stats]=regress(y,X);
            
            r2(ri,:,s)=stats(1);
            F(ri,:,s)=stats(2);
            p(ri,:,s)=stats(3);
            
            Y=X*b(ri,:,s)';
            Y=reshape(Y,roi_voxn,length(keptT));
            coupling(ri,:,s)=zeros(1,length(keptT));
            coupling(ri,:,s)=corr_col(gdata(:,keptT,s),Y);
            
            clear X
        end
        
    end
    couplingz=0.5*log((1+coupling)./(1-coupling));
    save([expdir '/' exp '/fmri/pattern_regression/' timeUnit '/roi/' froidir '/regression_LL2_lag' num2str(min(lags)) '-' num2str(max(lags))],'b','F','r2','p','couplingz','lags','rnames','keptT');
    clear b F p r2 rnames coupling
end

toc
