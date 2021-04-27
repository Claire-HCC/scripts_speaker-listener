clear all
set_parameters;

modelName='model_linear_noVolterra';

for ei=1:2;
    exp=exp_parameters.experiments{ei};
    voln=exp_parameters.voln_listener{ei};
    
    tb=readtable([expdir exp '/sound/transcripts_'  exp '_Claire.xlsx']);
    tb(strcmp(tb.word,'um'),:)=[];
    tb(strcmp(tb.word,'oh'),:)=[];
    tb(strcmp(tb.word,'ya'),:)=[];
    tb(strcmp(tb.word,'yeah'),:)=[];
    tb(tb.sentencei==0,:)=[];
    tb(tb.tmin_tr>voln,:)=[];
    
    names={'syllabel','word','sentence','paragraph','noPrLabel'};
    
    onsets=cell(1,length(names));
    offsets=cell(1,length(names));
    durations=cell(1,length(names));
    
    onsets{5}=min(tb.tmin(isnan(tb.paragraphi)));
    durations{5}=max(tb.tmax(isnan(tb.paragraphi)))-min(tb.tmin(isnan(tb.paragraphi)));
    
    tb(isnan(tb.paragraphi),:)=[];
    
    onsets{1}=[];
    durations{1}=[];
    for wi=1:length(tb.tmin);
        dur_syl=(tb.tmax(wi)-tb.tmin(wi))/tb.syllableN(wi);
        onsets{1}(end+1)=tb.tmin(wi);
        durations{1}(end+1)=dur_syl;
        
        if tb.syllableN(wi)>1;
            for syli=2:tb.syllableN(wi);
                onsets{1}(end+1)=tb.tmin(wi)+dur_syl*(syli-1);
                durations{1}(end+1)=dur_syl;
            end
        end
    end
    
    onsets{2}=tb.tmin;
    onsets{3}=tb.tmin(find([1; abs(diff(tb.sentencei))]>0  & tb.sentencei>0));
    onsets{4}=tb.tmin(find(([1; abs(diff(tb.paragraphi))]>0 |isnan([1; abs(diff(tb.paragraphi))]) & tb.paragraphi>0)));
    
    offsets{2}=tb.tmax;
    offsets{3}=tb.tmax(find([abs(diff(tb.sentencei)); 1]>0 & tb.sentencei>0));
    offsets{4}=tb.tmax(find(([abs(diff(tb.paragraphi)); 1]>0  | isnan([abs(diff(tb.paragraphi)); 1])& tb.paragraphi>0)));
    
    for levi=2:4;
        durations{levi}=offsets{levi}-onsets{levi};
    end
    
    
    pmod(1).name{1}='wd_sylN';
    pmod(2).name{1}='sn_wdN';
    pmod(3).name{1}='pr_snN';
    
    for lev=1:3;
        temp=zeros(1,length(onsets{lev}));
        for ui=1:length(onsets{lev+1});
            idx=find(onsets{lev}>=onsets{lev+1}(ui) & onsets{lev}<offsets{lev+1}(ui));
            temp(idx)=1:length(idx);
        end
      %  pmod(lev).param{1} = log(temp+1)-nanmean(log(temp+1));
       pmod(lev).param{1} = temp-nanmean(temp);
        pmod(lev).poly{1}  = 1;
    end
    
    mkdir([expdir exp '\fmri\temporal\stimulus2brain_glm\' modelName '\L1_spm\']);
    save([expdir exp '\fmri\temporal\stimulus2brain_glm\' modelName '\L1_spm\conditionSpec.mat'],'onsets','durations','names','pmod');
    
end


