NRFs={'linear', 'log', 'triangle','boxcar','linearR','expn','time'};
fsize=[16 4];
figure('unit','centimeter','position',[0 0 fsize]);
spi=1;
for nrf=1:3;%:length(NRFs);
    accumulationF=NRFs{nrf};
    
    temp=ones(1,11);
    if strcmp(NRF,'linear');
        temp=cumsum(temp);
    elseif     strcmp(NRF,'log');
        temp=cumsum(temp);
        temp=log(temp+1);
        
    elseif strcmp(NRF,'triangle');
        temp=cumsum(temp);
        mx=quantile(unique(temp),0.5);
        temp= mx-(abs(temp-mx));
        
    elseif strcmp(NRF,'boxcar');
        
    elseif strcmp(NRF,'linearR');
        temp=cumsum(temp);
        temp=fliplr(temp);
        
    elseif  strcmp(NRF,'expn');
        temp=cumsum(temp);
        temp=exp(temp+1);
        
    elseif strcmp(NRF,'time');
        temp=1:length(temp);
    end
    
    subplot(1,3,spi); spi=spi+1;
    xlim([0 17]);
    if strcmp(NRF,'boxcar');
        temp=[repmat(0,1,2) temp repmat(0,1,2)];
    else
        temp=[repmat(min(temp),1,2) temp repmat(min(temp),1,2)];
    end
    plot([1 2 2 3:12 12 13],temp,'color','k','linewidth',2);
    ylim([min(temp) max(temp)]);
    xlim([1 13])
    hold on
    line([0 length(temp)],[min(temp) min(temp)],'color','k','linewidth',2)
    axis off
    hold off
end