
clear all;

tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;

froidir='mor';

load([expdir '/roi_mask/' froidir '/roi_id_region.mat'],'roi_table');
tic % 15 min

% cropt start because ther eis clearly a spech-start effect in the
% lsiteners' data
crop_start=10;

lags=-4:4;
for ei=3;
    exp=experiments{ei};
    rnames=dir([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/'  froidir '/listenerAll_*.mat']);
    rnames={rnames.name};
    rnames=strrep(rnames,'listenerAll_','');
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
        
        load([expdir '/' exp '/fmri/timeseries/' timeUnit '/roi/' froidir '/listenerAll_' rname '.mat' ],'gdata');
        
          %% isc: average across voxels; 
        % then zscore overtime, so the results won't be dominated by
        % subjects with higher absolute values.
        gdata=zscore(nanmean(gdata,1),0,2);
        
        roi_voxn=size(gdata,1);
        roi_tn=size(gdata,2);
        keptT=(abs(min(lags))+1+crop_start):(roi_tn-max([lags,0]));
        
        for s=1:size(gdata,3);
            
            othersi=1:size(gdata,3);
            others=othersi(othersi~=s);
            others=nanmean(gdata(:,:,othersi),3);
            
            self=gdata(:,:,s);

            y=self;
            y=y(:,keptT);
            y=y(:);
            
            for li=1:length(lags);
                X(:,:,li)=others(:,keptT+lags(li));
            end
            
            X=reshape(X,roi_voxn*length(keptT),length(lags));
            
            % centralized X
            X=X-mean(X);
           
            % add an constant
            X=[ones(size(X,1),1) X];
            
            [b(ri,:,s),~,r(ri,:,s),~,stats]=regress(y,X);
            
            r2(ri,:,s)=stats(1);
            F(ri,:,s)=stats(2);
            p(ri,:,s)=stats(3);
            
             Y(ri,:)=X*b(ri,:,s)';
%             Y=reshape(Y,roi_voxn,length(keptT));
%             coupling(ri,:,s)=zeros(1,roi_tn);
%             coupling(ri,keptT,s)=corr_col(others(:,keptT),Y);
            
            clear X
        end
        
    end
  %  couplingz=0.5*log((1+coupling)./(1-coupling));
  b_ll=b;
  F_ll=F;
  r2_ll=r2;
  p_ll=p;
  r_ll=r;
  lags_ll=lags;
  rnames_ll=rnames;
  keptT_ll=keptT;
  Y_ll=Y;
  
    save([expdir '/' exp '/fmri/isc_regression/' timeUnit '/roi/' froidir '/regression_LL_lag' num2str(min(lags)) '-' num2str(max(lags))],'b_ll','F_ll','r2_ll','p_ll','r_ll','lags_ll','rnames_ll','keptT_ll','Y_ll');
    clear b F p r2 rnames coupling r Y
end

toc
