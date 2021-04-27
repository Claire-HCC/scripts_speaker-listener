% pause makes a huge difference, aotumatic detection of unvoiced pause and
% see whether story with more long pause tends to show clearere peaklag
% gradient. test with the ABC data.

% word rate =the output

% fit data
% slope after big break

% pallier vs. hierarchical
% mean higher activation in higher region vs. not realy


% play with the compression R distribution
% distribution change with level?

% mutual infor,ation can detect hierarchical structure?
% minus closed consitutent

% slower speech rate, longer lag
% log

% hrf convolution = low-pass filter?
% r stats
% transformation, irregular
close all
clear all
set_parameters;



crop_start=round(25*1.5);
crop_end=20*1.5;

exp='sherlock';
tb=readtable([expdir exp '/sound/transcripts_'  exp '_Claire.xlsx']);

lags=(-20*1.5):(20*1.5);
permN=1000;
simN=100;

levn=6;
accumulationFs={'linear', 'log', 'square'};
pause=1;
wdN=3000;
spi=1;
for u=3;%2:4
    for vi=0.4;%[0.1:0.3:1.5];
        for speechR=1;%[0.5 1 2]; % 1 is the speech rate of sherlock;
            for af=1;%1:3;
                vs_hrf_sec=[];
                for simi=1%:simN;
                    
                    accumulationF=accumulationFs{af};
                    dur_wd=round((tb.tmax-tb.tmin)*1000/speechR);
                    
                    mu = log((u^2)/sqrt(vi+u^2));
                    sigma = sqrt(log(vi/(u^2)+1));
                    compressionR=round(lognrnd(mu,sigma,1000,1));
                    compressionR(compressionR<0)=[];
                    
                    % create word onset vectors
                    dur_wd_samples=[];
                    dur_syl_samples=[];
                    while length(dur_wd_samples)<wdN;
                        samplei=datasample(1:length(dur_wd),1);
                        dur_wd_samples(end+1)=dur_wd(samplei);
                        sylN=tb.syllableN(samplei);
                        
                        % make sure that syllable duration is intager. Put
                        % additional duration in the 1st syllable.
                        dur_otherSyl=floor(dur_wd(samplei)/sylN);
                        dur_1stSyl=dur_wd(samplei)-(dur_otherSyl*(sylN-1));
                        
                        dur_syl_samples(end+1)=dur_1stSyl;
                        dur_syl_samples((end+1):(end+sylN-1))=dur_otherSyl;
                    end
                    
                    tn=sum(dur_wd_samples);
                    onsetsVectors=zeros(tn,levn+1);
                    onsetsVectors(cumsum(round(dur_syl_samples))+1,1)=1;
                    onsetsVectors(cumsum(dur_wd_samples)+1,2)=1;
                    
                    for lev=3:(levn+1);
                        n_before=sum(onsetsVectors(:,lev-1));
                        onsets_before=find(onsetsVectors(:,lev-1)==1);
                        compressionR_sample=datasample(compressionR,wdN);
                        n_current=sum(cumsum(compressionR_sample)<n_before);
                        compressionR_sample=compressionR_sample(1:n_current);
                        onsets_current=onsets_before(cumsum(compressionR_sample));
                        onsetsVectors(onsets_current,lev)=1;
                    end
                    
                    onsetsVectors(1,:)=1;
                    onsetsVectors=onsetsVectors(1:tn,:);
                    
                    % remove levels with less than 2 units
                    onsetsVectors(:,sum(onsetsVectors)<2)=[];
                    % syllable level does not count
                    levn_temp=size(onsetsVectors,2);
                    
                    v=[];
                    for lev=2:(levn_temp);
                        onsets_current=find(onsetsVectors(:,lev)==1);
                        
                        for onseti=1:length(onsets_current);
                            if onseti==length(onsets_current);
                                idx=onsets_current(onseti):tn;
                            else
                                idx=onsets_current(onseti):(onsets_current(onseti+1)-1);
                            end
                            temp=onsetsVectors(idx,lev-1);
                            
                            % accumulation function
                            % the size of pause effect might make a difference
                            if strcmp(accumulationF,'linear');
                                temp=cumsum(temp);
                            elseif     strcmp(accumulationF,'log');
                                temp=cumsum(temp);
                                temp=log(temp+1);
                            elseif strcmp(accumulationF,'square');
                                temp=cumsum(temp);
                                temp=(max((temp-mean(temp))+1).^2)-(temp-mean(temp)).^2;
                            end
                            
                            v(idx,lev)=temp;
                        end
                    end
                    
                    
                    % insert pause
                    if pause==1;
                        onsets=find(onsetsVectors(:,levn_temp)==1);
                        v_nopause=v;
                        v=[];
                        for ci=1:(length(onsets));
                            if ci==length(onsets);
                                idx=onsets(ci):tn;
                            else
                                idx=onsets(ci):(onsets(ci+1)-1);
                            end
                            v=[v; zeros(3*1000,levn_temp) ;v_nopause(idx,:)];
                        end
                    end
                    
                    % hrf convolution
                    xBF.dt=0.001;
                    xBF.name='hrf (with time and dispersion derivatives)';
                    bf = spm_get_bf(xBF);
                    
                    v_hrf=[];
                    for xi=1:size(v,2);
                        temp=conv(v(:,xi),bf.bf(:,1));
                        v_hrf(:,xi)=temp;
                    end
                    
                    % downsampling to  1 sec resolution
                    v_hrf_sec=resample(v_hrf,1,1000);
                    v_hrf_sec=v_hrf_sec((crop_start+1):(end-crop_end),:);
                    % remove syllable level
                    v_hrf_sec(:,1)=[];
                    
                    vs_hrf_sec(:,:,simi)=v_hrf_sec;
                    
                    levn_temp=size(v_hrf_sec,2);
                    
                    % compute lagcorrelation
                    r=[];
                    tn=size(v_hrf_sec,1);
                    t_real=round((tn-1)/2)+1;
                    for sdi=1:levn_temp;
                        for tgi=1:levn_temp;
                            r(sdi,tgi,:)=circularlagcorr(v_hrf_sec(:,sdi),v_hrf_sec(:,tgi),[-(t_real-1):(t_real-1)]);
                        end
                    end
                    
                    % stats
                    % -0.00001 so to avoid inf after atanh
                    rz=atanh(r-0.000001);
                    
                    ts_shift=1:tn;
                    ts_shift=ts_shift((ts_shift+min(lags))>=1 & (ts_shift+max(lags))<=tn);
                    peaks_shift=[];
                    for perm=1:permN;
                        ti=randi(length(ts_shift));
                        t_shift=ts_shift(ti);
                        [peaks_shift(:,:,perm),lagi]=max(rz(:,:,t_shift+lags),[],3);
                    end
                    rz=rz(:,:,t_real+lags);
                    p=mean(permute(repmat(peaks_shift,1,1,1,length(lags)),[1 2 4 3])>rz,4);
                    pfwe=p*length(p(:))/length(lags);
                    
                    peaks=nan(levn,levn);
                    peakLags=nan(levn,levn);
                    [peaks(1:levn_temp,1:levn_temp),lagi]=max(rz,[],3);
                    peakLags(1:levn_temp,1:levn_temp)=lags(lagi);
                    peakLags_pfwe=peakLags;
                    for sdi=1:levn_temp;
                        for tgi=1:levn_temp;
                            if squeeze(pfwe(sdi,tgi,lagi(sdi,tgi)))>.05;
                                peakLags_pfwe(sdi,tgi)=NaN;
                            end
                        end
                    end
                    peaks_pfwe=peaks;
                    peaks_pfwe(isnan(peakLags_pfwe))=NaN;
                    
                    % remove syllable level
                    v(:,1)=[];
                    onsetsVectors(:,1)=[];
                    
                    
                end
            end
            %
        end
    end
end

% 
% 
% fsize=[27 17];
% set(gcf,'unit','centimeter','position',[0 0 fsize]);
% % peak lag imagesc
% 
% subplot(3,5,spi);
% 
% %            peakLags_bined=discretize(peakLags,[-20 -15 -10 -6 -3 0 1 4  7  11 16 20]*1.5);
% %    im=imagesc(peakLags_bined,[2 10]);
% 
% im=imagesc(peakLags,[-5 5]);
% colormap jet
% set(gca,'color','k')
% set(im,'AlphaData',~isnan(peakLags_pfwe));
% wpm=round(length(dur_wd_samples)/(sum(dur_wd_samples)/60000));
% 
% title({'Compression ratio:',sprintf('Mean=%d, SD=%.2f',u,sqrt(vi))});
% % title({sprintf('%d word per min',wpm)});
% 
% spi=spi+1;

cols=jet(7);

figure
for lev=1:levn;
    subplot(2,1,1);
    plot(v(:,lev),'color',cols(lev+1,:));
    hold on;
    xlim([0 size(v,1)])
    set(gca,'xticklabel',[])
    ylabel('Neural activity')
    set(gca,'fontsize',14)
    subplot(2,1,2);
    plot(v_hrf_sec(:,lev),'color',cols(lev+1,:));
    hold on;
    xlim([0 size(v_hrf_sec,1)])
    ylabel('BOLD signal')
    set(gca,'fontsize',14)
end
xlabel('Time (sec)');
legend(cellstr(num2str([1:6]')),'orientation','horizontal')
legend boxoff

% 
% figure
% for sdi=1:levn;
%     subplot(2,3,sdi);
%     imagesc(squeeze(zscore(r(sdi,:,:),0,3)));
% end
% 
% 
% figure
% for sdi=1:levn;
%     subplot(2,3,sdi);
%     for tgi=1:levn;
%         
%         plot(lags,squeeze(rz(sdi,tgi,:)),'color',cols(tgi+1,:));
%         hold on
%     end
%     hold off
% end
% 
% figure('unit','centimeter','position',[0 0 11 12]);
% histogram(compressionR,length(unique(compressionR)),'normalization','probability')
% xlabel({'Compression ratio','(e.g. word per sentence)'})
% ylabel('Probability (%)')
% set(gca,'fontsize',16)
% title({'Lognormal distribution',sprintf('Mean=%d, SD=%.2f',u,sqrt(vi))});

% frequency analysis
Sxx=[];
Syy=[];
Sxy=[];
for sdi=1:6;
    for tgi=1:6;
        [Sxx(sdi,tgi,:) freq] = pwelch(zscore(v_hrf_sec(:,sdi)),100,50,[],1);
        [Syy(sdi,tgi,:) freq] = pwelch(zscore(v_hrf_sec(:,tgi)),100,50,[],1);
        [Sxy(sdi,tgi,:) freq] =cpsd(zscore(v_hrf_sec(:,sdi)),zscore(v_hrf_sec(:,tgi)),100,50,[],1);

    end
end
coherences=abs(Sxy).^2./(Sxx.*Syy);
phases=angle(Sxy);

cols=jet(7);

%¡@psd
figure;
for ni=1:6;
    plot((freq(2:end)),squeeze(Sxx(ni,ni,2:end)),'color',cols(ni+1,:),'linewidth',2);
    hold on
end
hold off
xlabel('Frequency (Hz)');
ylabel('PSD (1/Hz)');
xlim(log([freq(2) freq(end)]))
set(gca,'xtick',([0.01 0.04 0.1 0.33]));
set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
grid on
set(gca,'fontsize',14);
 title({'Compression ratio:',sprintf('Mean=%d, SD=%.2f',u,sqrt(vi))});
legend(cellstr(num2str([1:6]')))
legend boxoff

% cross-network coherence
figure;
for tgi=1:6;
    plot((freq(2:end)),squeeze(coherences(6,tgi,2:end)),'color',cols(tgi+1,:),'linewidth',2);
    hold on
end
hold off
xlabel('Frequency (Hz)');
ylabel(' cross-network Coherence (1/Hz)');
xlim(([freq(2) freq(end)]))
set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
grid on
set(gca,'fontsize',14)
 title({'Compression ratio:',sprintf('Mean=%d, SD=%.2f',u,sqrt(vi))});

% cross-network phase
figure;
for tgi=1:6;
    plot((freq(2:end)),squeeze(phases(6,tgi,2:end)),'color',cols(tgi+1,:),'linewidth',2);
    hold on
end
hold off
xlabel('Frequency (Hz)');
ylabel('Cross-network phase');
xlim(log([freq(2) freq(end)]))
set(gca,'xtick',log([0.01 0.04 0.1 0.33]));
set(gca,'xticklabels',[0.01 0.04 0.1 0.33]);
grid on
set(gca,'fontsize',14)
 title({'Compression ratio:',sprintf('Mean=%d, SD=%.2f',u,sqrt(vi))});

