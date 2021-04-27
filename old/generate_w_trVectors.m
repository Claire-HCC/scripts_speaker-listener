clear all
set_parameters;
exp='merlin';
w_table=readtable([expdir exp '/sound/wordTimeStamp.csv']);

ws={{'na'},{'na','um'}};
trVector=table();
for wi=1:length(ws);
    w=ws{wi};
    
    i=ismember(w_table.word,w);;
    temp=diff(i);
    
    s=find(temp==1)+1;
    e=find(temp==-1);
    if length(s)<length(e);
        s=[1; s];
    elseif length(s)>length(e);
        e=[1; e];
    end
    
    timeStamp=([w_table.trmin(s) w_table.trmax(e)]);
    
    % only keep keywords that last longer or equal to  1 tr
    timeStamp((timeStamp(:,2)-timeStamp(:,1))<1,:)=[];
    
    timeStamp=round(timeStamp);
    
    % starting with tr 1
    timeStamp(timeStamp==0)=1;
    
    w_trVectors(1:round(max(w_table.trmax)),strjoin(w,'_'))=table(0);
    for ti=1:size(timeStamp,1);
        w_trVectors(timeStamp(ti,1):timeStamp(ti,2),strjoin(w,'_'))=table(1);
    end
end

save([expdir exp '/sound/w_trVectors.mat'],'w_trVectors')