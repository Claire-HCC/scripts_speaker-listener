clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-3;
binSize=1;

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};

for ei=3%1:4;
    exp=experiments{ei};
    
    %   mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/textVectors/']);
   
     load([expdir exp '/sound/text_trVectors_binSize' num2str(binSize)   '.mat'],'text_trVectors');
    textvarNames=text_trVectors.Properties.VariableNames;
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
    r2_sl=r2;
    r2_sl_m=nanmean(r2_sl,3);
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags))],'r2','rnames');
    r2_ll=r2;
    r2_ll_m=nanmean(r2,3);
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/LL_minus1Leave1out_bined/binSize' num2str(binSize) '_lag0-0.mat'],'r2');
    r2_ll0=nanmean(r2,4);
    r2_ll0_m=nanmean(r2_ll0,3);
    
    keptT=[min(find(~isnan(r2_sl(1,:,1)))):max(find(~isnan(r2_sl(1,:,1))))];
    [~,tn,listenerN]=size(r2_sl);
    
    % to get the right degree of freedom
  %   keptT=min(keptT):binSize:max(keptT);
    
    r2=r2_ll0_m;
    ris=find(sum(isnan(r2(:,keptT)),2)==0);
   r=nan([length(rnames) 31]);
   fsize=[20 20];
    figure('unit','centimeter','position',[0 0 fsize]);
    for vi=1:length(textvarNames);
        textvarName=textvarNames{vi};
        r(ris,:)=lagcorr_claire(repmat(text_trVectors.(textvarName)(keptT),1,length(ris)),r2(ris,keptT)',-15:15)';
        subplot(3,3,vi);
        %        % imagesc(r',[-0.1 0.1]);
           %   plot(-10:10,r);
        
        rsig=find(sum(abs(r)>0.15,2)>0);
        if ~isempty(rsig);
        plot(-15:15,r(rsig,:),'linewidth',2);
        legend(rnames(rsig));
        legend boxoff
        end
         grid on
        title(textvarName);
    end
    
    %     save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/textVectors/LL_minus1Leave1out_binSize' num2str(binSize) '_lag0-0.mat'],...
    %         'r_eventb','p_eventb','pfdr_eventb','r_sentenceb','p_sentenceb','pfdr_sentenceb','r_wrate','p_wrate','pfdr_wrate','r_naum','p_naum','pfdr_naum');
    %
    %
    
    %     save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/textVectors/SLvsLL_leave1out_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],...
    %         'r_eventb','p_eventb','pfdr_eventb','r_sentenceb','p_sentenceb','pfdr_sentenceb','r_wrate','p_wrate','pfdr_wrate','r_naum','p_naum','pfdr_naum');
    %
end


