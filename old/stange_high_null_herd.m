for i =1:1000;
    %% speaker has its unique component and that is also in the listeners data
    %     speaker_component=rand(100,500);
    %     listener_component=rand(100,500);
    %     listeners=rand(100,500,18)+0.1*repmat(speaker_component,1,1,18)+repmat(listener_component,1,1,18);
    % speaker=speaker_component+rand(100,500);
    
    % listener is not different from the listeners
    shared_component=rand(100,500);
    listener=rand(100,500);
    listeners=repmat(shared_component,1,1,18)+rand(100,500,18);
    speaker=0.5*shared_component+rand(100,500);
    
    sl(:,1)=corr_col(speaker,mean(listeners,3));
    for s=1:18;
        othersi=1:18;
        othersi=othersi(othersi~=s);
        ll(:,s)=corr_col(listeners(:,:,s),mean(listeners(:,:,othersi),3));
    end
    h(i)=corr(real(atanh(sl)),mean(real(atanh(ll)),2));
end


% sl ll caculated with two different groups of subjects
clear h
for i =1:1000;
    shared_component=rand(100,500);
    listener=rand(100,500);
    listeners=repmat(shared_component,1,1,18)+rand(100,500,18);
    listenres=zscore(listeners,0,2);
    speaker=0.5*shared_component+rand(100,500);
    speaker=zscore(speaker,0,1);
    
    sl(:,1)=corr_col(speaker,mean(listeners(:,:,1:9),3));
    
    for s=10:18;
        self=listeners(:,:,s);
        othersi=10:18;
        othersi=othersi(othersi~=s);
        others=nanmean(listeners(:,:,othersi),3);
        ll(:,s)=corr_col(self,others);
    end
    h(i)=corr(real(atanh(sl)),mean(real(atanh(ll)),2));
end

