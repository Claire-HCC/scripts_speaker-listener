
clear all;

tic
loc='cluster';
set_parameters;
timeUnit='tr' ;

lags=-5:5;
for ei=1%:2;
    exp=experiments{ei};
    b_SL=[];
    
    load([expdir '/' exp '/fmri/timeseries/' timeUnit '/wholeBrain//zscore_listenerAll_' rname ],'gdata','keptvox');
    
    for vi=1:length(keptvox);
        for s=1:size(gdata,3);
            
            self=gdata(vi,:,s);
            othersi=1:size(gdata,3);
            othersi(othersi==s)=[];
            others=nanmean(gdata(vi,:,othersi),3);
            
            keptT=(max(lags)+1):(size(data,2)-max(lags));
            self=self(:,keptT);
            
            for li=1:length(lags);
                data_mat(:,:,li)=data(vi,keptT+lags(li));
            end
            
            self=self(:);
            
            data_mat=reshape(data_mat,size(data_mat,1)*size(data_mat,2),length(lags));
            b_SL(vi,:,s)=ridge(self,data_mat,0);
            clear data_mat
        end
        
    end
end

save([expdir '/' exp '/fmri/pattern_ridge/' timeUnit '/wholeBrain/beta_SL' ],'b_SL','lags');
% clear b_SL roi_ids rnames

toc
