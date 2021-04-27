clear all
close all
set_parameters;
timeUnit='tr' ;
froidir='mor';
binSize=10;
lags=-10:-3;
load([expdir '/roi_mask/mor/' 'roi_id_region.mat'],'roi_table');
exp=experiments{2};

for ei=3%1:4;
    exp=experiments{ei};
    
    load([expdir '/' exp '/fmri/pattern/regression/' timeUnit '/roi/' froidir '/SL_leave1out_bined/binSize' num2str(binSize) '_lag' num2str(min(lags)) '-' num2str(max(lags)) ],'r2','rnames');
    r2_sl=r2;
    
    keptT=[min(find(~isnan(r2_sl(1,:,1)))):max(find(~isnan(r2_sl(1,:,1))))];
    [~,tn,listenerN]=size(r2_sl);
    
    w_table=readtable([expdir exp '/sound/wordTimeStamp.csv']);
    load([expdir exp '/sound/w_trVectors.mat'],'w_trVectors');
    w_trVectors=w_trVectors(1:tn,:);
    
    na=zeros(tn,1);
    naum=zeros(tn,1);
    wrate=zeros(tn,1);
    speech=zeros(tn,1);
    visual=zeros(tn,1);
    sentenceb=zeros(tn,1);
    mainCharacterb=zeros(tn,1);
    arousal=zeros(tn,1);
    sentenced=[0 ;diff(w_table.sentence)];
    eventb=zeros(tn,1);
    eventd=[ 0 ;diff(w_table.event)];
    mainCharacterd=[ 0 ;diff(w_table.mainCharacter)];
    
    for t=1:tn;
        t_bin=t:(t+binSize-1);
        
        if min(t_bin)>=1 & max(t_bin)<=tn;
            wi=(ismember(round(w_table.trmin),t_bin));% & ismember(round(w_table.trmax),t_bin) );
            na(t)=sum(table2array(w_trVectors(t_bin,'na')));
            naum(t)=sum(table2array(w_trVectors(t_bin,'na_um')));
            wrate(t)=sum(wi & ~ismember(w_table.word,{'na','um'}));
            sentenceb(t)=sum(sentenced(wi));
            eventb(t)=sum(eventd(wi));
            speech(t)=sum(w_table.speech(wi));
            visual(t)=sum(w_table.visual(wi));
            arousal(t)=sum(w_table.arousal(wi));
            mainCharacterb(t)=sum(mainCharacterd(wi));
        end
    end
    text_trVectors=table(eventb,sentenceb,wrate,na,naum,visual,speech,mainCharacterb,arousal);
    save([expdir exp '/sound/text_trVectors_binSize' num2str(binSize)   '.mat'],'text_trVectors');
end



