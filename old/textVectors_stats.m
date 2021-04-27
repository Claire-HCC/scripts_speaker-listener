clear all;
close all
tic
% loc='cluster';
set_parameters;
timeUnit='tr' ;
froidir='mor';
lags=-10:-3;
binSize=10;

load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};

for ei=3%1:4;
    
    exp=experiments{ei};
    mkdir([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/textVectors/']);
    load([expdir exp '/sound/text_trVectors_binSize' num2str(binSize)   '.mat'],'text_trVectors');
    
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
    keptT=min(keptT):binSize:max(keptT);
    
    
    r2=r2_ll0_m;
    r=nan([length(rnames) size(text_trVectors,2)]);
    p=nan([length(rnames) size(text_trVectors,2)]);
    pfdr=nan([length(rnames) size(text_trVectors,2)]);
    textvarNames=text_trVectors.Properties.VariableNames;
    for vi=1:length(textvarNames);
        textvarName=textvarNames{vi};
        for ri=1:length(rnames);
            [r(ri,vi), p(ri,vi)]=corr(text_trVectors.(textvarName)(keptT),r2(ri,keptT)');
        end
        [~,~,pfdr(~isnan(r(:,vi)),vi)]=fdr(p(~isnan(p(:,vi)),vi));
    end
    
    %     save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/textVectors/LL_minus1Leave1out_binSize' num2str(binSize) '_lag0-0.mat'],...
    %         'r_eventb','p_eventb','pfdr_eventb','r_sentenceb','p_sentenceb','pfdr_sentenceb','r_wrate','p_wrate','pfdr_wrate','r_naum','p_naum','pfdr_naum');
    %
    
    %     r2=r2_sl_m-r2_ll_m;
    %     r=nan([length(rnames) size(text_trVectors,2)]);
    %     p=nan([length(rnames) size(text_trVectors,2)]);
    %     pfdr=nan([length(rnames) size(text_trVectors,2)]);
    %     for vi=1:size(text_trVectors,2);
    %         for ri=1:length(rnames);
    %             [r(ri,vi), p(ri,vi)]=corr(text_trVectors(keptT,vi),r2_ll0_m(ri,keptT)');
    %         end
    %         [~,~,pfdr(~isnan(r(:,vi)))]=fdr(p_sentenceb(~isnan(p(:,vi))));
    %     end
    
    %     save([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/textVectors/SLvsLL_leave1out_binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) '.mat'],...
    %         'r_eventb','p_eventb','pfdr_eventb','r_sentenceb','p_sentenceb','pfdr_sentenceb','r_wrate','p_wrate','pfdr_wrate','r_naum','p_naum','pfdr_naum');
    %
end


