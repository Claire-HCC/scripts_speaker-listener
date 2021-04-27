

close all
clear all
loc='mypc';
set_parameters;

crop_start=round(25);
crop_end=20;

exp='sherlock';
tb=readtable([expdir exp '/sound/transcripts_'  exp '_Claire.xlsx']);

% get hrf
 xBF.dt=0.001; 
 xBF.name='hrf (with time and dispersion derivatives)';
bf = spm_get_bf(xBF);
% load([expdir '/simulation/bf.mat']);

permN=1000;
lags=-20:20;
simN=1;
levn=6;
NRFs={'linear', 'log', 'square'};
pauses_sec=[0  1.5 3];
speechRs=[0.5 1 2];
wdN=500;
spi=1;

for u=3%:4
    for vi=0;%0.4;%[0.1:0.3:1.5];
        mu = log((u^2)/sqrt(vi+u^2));
        sigma = sqrt(log(vi/(u^2)+1));
        constituentL=round(lognrnd(mu,sigma,1000,1));
        constituentL(constituentL<0)=[];
        constituentL_sample=datasample(constituentL,wdN);
        
        for psi=3;%1:length(pauses_sec);
            pause_sec=pauses_sec(psi);
            
            for speechRi= 1%:length(speechRs); % 1 is the speech rate of sherlock;
                speechR=speechRs(speechRi);
                dur_wd=round(285/speechR);
                
                for af=1%:3;
                    NRF=NRFs{af};
                    % create word onset vectors
                    dur_wd_samples=repmat(dur_wd,wdN,1);
                    dur_syl_samples=dur_wd_samples;
                    
                    tn=sum(dur_wd_samples);
                    onsetsVectors=zeros(levn+1,tn);
                    onsetsVectors(1,cumsum(round(dur_syl_samples))+1)=1;
                    onsetsVectors(2,cumsum(dur_wd_samples)+1)=1;
                    
                    for lev=3:(levn+1);
                        n_below=sum(onsetsVectors(lev-1,:));
                        onsets_below=find(onsetsVectors(lev-1,:)==1);
                        
                        n_current=sum(cumsum(constituentL_sample)<=n_below);
                        constituentL_sample_temp=constituentL_sample(1:n_current);
                        onsets_current=onsets_below(cumsum(constituentL_sample_temp));
                        onsetsVectors(lev,onsets_current)=1;
                    end
                    
                    onsetsVectors(:,1)=1;
                    onsetsVectors=onsetsVectors(:,1:tn);
                    
                    % remove levels with less than 2 units
                    onsetsVectors(sum(onsetsVectors,2)<2,:)=[];
                    % syllable level does not count
                    levn_temp=size(onsetsVectors,1);
                    
                    v=[];
                    for lev=2:(levn_temp);
                        onsets_current=find(onsetsVectors(lev,:)==1);
                        
                        for onseti=1:length(onsets_current);
                            if onseti==length(onsets_current);
                                idx=onsets_current(onseti):tn;
                            else
                                idx=onsets_current(onseti):(onsets_current(onseti+1)-1);
                            end
                            temp=onsetsVectors(lev-1,idx);
                            
                            % accumulation function
                            % the size of pause effect might make a difference
                            if strcmp(NRF,'linear');
                                temp=cumsum(temp);
                            elseif     strcmp(NRF,'log');
                                temp=cumsum(temp);
                                temp=log(temp+1);
                            elseif strcmp(NRF,'square');
                                temp=cumsum(temp);
                                temp=(max((temp-mean(temp))+1).^2)-(temp-mean(temp)).^2;
                            end
                            
                            v(lev,idx)=temp;
                        end
                    end
                    
                    % insert pause
                    % the effect size of pause might make a difference
                    onsets=find(onsetsVectors(levn_temp,:)==1);
                    v_nopause=v;
                    v=[];
                    onsetsVectors_nopause=onsetsVectors;
                    onsetsVectors=[];
                    for ci=1:(length(onsets));
                        if ci==length(onsets);
                            idx=onsets(ci):tn;
                        else
                            idx=onsets(ci):(onsets(ci+1)-1);
                        end
                        v=[v zeros(levn_temp,pause_sec*1000) v_nopause(:,idx)];
                        onsetsVectors=[onsetsVectors zeros(levn_temp,pause_sec*1000) onsetsVectors_nopause(:,idx)];
                    end
                    
                    % hrf convolution
                    v_hrf=[];
                    for xi=1:size(v,1);
                        temp=conv(zscore(v(xi,:)),bf.bf(:,1));
                        v_hrf(xi,:)=temp;
                    end
                    
                    
                    % downsampling to  1 sec resolution
                    v_hrf_sec_temp=resample(v_hrf',1,1000)';
                    % remove syllable level
                    v_hrf_sec_temp(1,:)=[];
                    
                    v_hrf_sec_temp=v_hrf_sec_temp(:,1:floor(size(v,2)/1000));
                    v_hrf_sec_temp=v_hrf_sec_temp(:,(crop_start+1):end);
                    
                    v_hrf_sec=v_hrf_sec_temp;
                    
                    [levn_temp, tn]=size(v_hrf_sec);
                    
                    v= v(:,(crop_start*1000+1):(crop_start+size(v_hrf_sec_temp,2))*1000);
                    onsetsVectors=onsetsVectors(:,(crop_start*1000+1):(crop_start+size(v_hrf_sec_temp,2))*1000);
                    % remove syllable level
                    v(1,:)=[];
                    onsetsVectors(1,:)=[];
                    
                    % compute lagcorrelation
                    r=[];
                    tn=size(v_hrf_sec,2);
                    t_real=round((tn-1)/2)+1;
                    for sdi=1:levn_temp;
                        for tgi=1:levn_temp;
                            r(sdi,tgi,:)=circularlagcorr(v_hrf_sec(sdi,:),v_hrf_sec(tgi,:),[-(t_real-1):(t_real-1)]);
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
                    peakLags_pfwe
                    
                    %                     figure;
                    %                     im=imagesc(peakLags_pfwe,[-5 5]);
                    %                     colormap jet
                    %                     set(gca,'color','k')
                    %                     set(im,'AlphaData',~isnan(peakLags_pfwe));
                    %
                    cols=jet(7);
                    subplot(3,1,af);
                    for lev=1:levn;
                        %                         subplot(2,1,1);
                        %                         plot(v(lev,:),'color',cols(lev+1,:));
                        %                         hold on;
                        %                         xlim([0 size(v,2)])
                        %                         set(gca,'xticklabel',[])
                        %                         ylabel('Neural activity')
                        %                         set(gca,'fontsize',14)
                        %                         subplot(2,1,2);
                        plot(v_hrf_sec(lev,:),'color',cols(lev+1,:));
                        hold on;
                        xlim([0 size(v_hrf_sec,2)])
                        ylabel('BOLD signal')
                        set(gca,'fontsize',14)
                    end
                    xlabel('Time (sec)');
                    legend(cellstr(num2str([1:6]')),'orientation','horizontal')
                    legend boxoff
                    %                  title(sprintf('speech rate: %.1f',speechR))
                    
                end
            end
        end
    end
end



