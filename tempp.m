
us=2; % 2:4; 3000 words are often not enough to create level 6 with more than to 2 consituents.
vars=2:3;
speechRs=1;
accumulationFs={'linear'}; % {'linear', 'log','expn', 'triangle','time'};
pauseLens=3;


loc='mypc';
set_parameters;
experiment='sherlock';
tb=readtable([expdir experiment '/sound/transcripts_'  experiment '_Claire.xlsx']);
load([expdir '/simulation/bf.mat']);
tr=1.5;
simN=30;
levn=6;
wdN=3000;

for ui=1:length(us);
    u=us(ui);
    
    for vari=1:length(vars);
        var=vars(vari);
        
        for speechRi=1:length(speechRs);
            speechR=speechRs(speechRi);
            
            for af=1:length(accumulationFs);
                accumulationF=accumulationFs{af};
                
                for pauseLeni=1:length(pauseLens);
                    pauseLen=pauseLens(pauseLeni);
                    
                    for pauseEffecti=1:length(pauseEffects);
                        pauseEffect=pauseEffects(pauseEffecti);
                        
                        if exist(sprintf('%s/simulation/sim%03d/pause%.2f_pauseL%.1f_m%d_v%.2f_sr%.2f_%s.mat',expdir,simN,pauseEffect,pauseLen,u,var,speechR, accumulationF))==0;
                           
                                for simi=1:simN;
                                    dur_wd=round((tb.tmax-tb.tmin)*1000/speechR);
                                    
                                    % speech R also applies to boundary
                                    % duration
                                    dur_boundary_sample=round(normrnd(pauseLen*1000,(pauseLen*1000)/6,1000,1)/speechR);
                                    dur_boundary_sample(dur_boundary_sample<0)=[];
                                    
                                    % parameters of the associated normal distribution
                                    mu = log((u^2)/sqrt(var+u^2));
                                    sigma = sqrt(log(var/(u^2)+1));
                                    
                                    compressionR=round(lognrnd(mu,sigma,1000,1));
                                    compressionR(compressionR<=0)=[];
                                    
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
                                    onsetsVectors_temp=zeros(levn+1,tn);
                                    onsetsVectors_temp(1,cumsum(round(dur_syl_samples))+1)=1;
                                    onsetsVectors_temp(2,cumsum(dur_wd_samples)+1)=1;
                                    
                                    for lev=3:(levn+1);
                                        n_below=sum(onsetsVectors_temp(lev-1,:));
                                        onsets_below=find(onsetsVectors_temp(lev-1,:)==1);
                                        compressionR_sample=datasample(compressionR,wdN);
                                        n_current=sum(cumsum(compressionR_sample)<=n_below);
                                        compressionR_sample=compressionR_sample(1:n_current);
                                        onsets_current=onsets_below(cumsum(compressionR_sample));
                                        onsetsVectors_temp(lev,onsets_current)=1;
                                    end
                                    
                                    onsetsVectors_temp(:,1)=1;
                                    onsetsVectors_temp=onsetsVectors_temp(:,1:tn);
                                    
                                    % remove levels with less than 2 units
                                    onsetsVectors_temp(sum(onsetsVectors_temp,2)<2,:)=[];
                                    % syllable level does not count
                                    levn_temp=size(onsetsVectors_temp,1);
                                    
                                    v_temp=[];
                                    for lev=2:(levn_temp);
                                        onsets_current=find(onsetsVectors_temp(lev,:)==1);
                                        
                                        for onseti=1:length(onsets_current);
                                            if onseti==length(onsets_current);
                                                idx=onsets_current(onseti):tn;
                                            else
                                                idx=onsets_current(onseti):(onsets_current(onseti+1)-1);
                                            end
                                            temp=onsetsVectors_temp(lev-1,idx);
                                            
                                            % accumulation function
                                            % the size of pause effect might make a difference
                                            if strcmp(accumulationF,'linear');
                                                temp=cumsum(temp);
                                            elseif     strcmp(accumulationF,'log');
                                                temp=cumsum(temp);
                                                temp=log(temp+1);
                                            elseif  strcmp(accumulationF,'expn');
                                                temp=cumsum(temp);
                                                temp=exp(temp+1);
                                            elseif strcmp(accumulationF,'triangle');
                                                temp=cumsum(temp);
                                                mx=quantile(unique(temp),0.5);
                                                temp= mx-(abs(temp-mx));
                                            elseif strcmp(accumulationF,'time');
                                                temp=1:length(idx);
                                            end
                                            
                                            v_temp(lev,idx)=temp;
                                        end
                                    end
                                    
                                    % insert pause
                                    % the effect size of pause might make a difference
                                    onsets=find(onsetsVectors_temp(levn_temp,:)==1);
                                    v_nopause=v_temp;
                                    v_temp=[];
                                    onsetsVectors_nopause=onsetsVectors_temp;
                                    onsetsVectors_temp=[];
                                    samplei=datasample(1:length(dur_boundary_sample),length(onsets));
                                    pauses=dur_boundary_sample(samplei);
                                    for ci=1:(length(onsets));
                                        if ci==length(onsets);
                                            idx=onsets(ci):tn;
                                        else
                                            idx=onsets(ci):(onsets(ci+1)-1);
                                        end
                                        mn=min(v_nopause,[],2);
                                        sd=std(v_nopause,[],2);
                                        pause_temp=repmat(mn-sd*pauseEffect,1,pauses(ci));
                                        v_temp=[v_temp pause_temp v_nopause(:,idx)];
                                        onsetsVectors_temp=[onsetsVectors_temp zeros(size(pause_temp)) onsetsVectors_nopause(:,idx)];
                                    end
                                    
                                    % remove syllable level
                                    v_temp(1,:)=[];
                                    onsetsVectors_temp(1,:)=[];
                                    
                                    % hrf convolution
                                    v_hrf_temp=[];
                                    for xi=1:size(v_temp,1);
                                        temp=conv(v_temp(xi,:),bf.bf(:,1));
                                        v_hrf_temp(xi,:)=temp;
                                    end
                                    
                                    % downsampling to  1 tr resolution
                                    v_hrf_tr_temp=v_hrf_temp(:,[0:(1000*tr):(end-tr*1000)]+round((tr*1000)/2));
                                    v_tr_temp=v_temp(:,[0:(1000*tr):(end-tr*1000)]+round((tr*1000)/2));
                                    v_temp=v_temp(:,[0:(end-tr*1000)]+round((tr*1000)/2));
                                    v_hrf_temp=v_hrf_temp(:,[0:(end-tr*1000)]+round((tr*1000)/2));
                                    onsetsVectors_temp=onsetsVectors_temp(:,[0:(end-tr*1000)]+round((tr*1000)/2));
                                    %       v_hrf_tr_temp=resample(v_hrf',1,tr*1000)';
                                    %      v_hrf_tr_temp=v_hrf_tr_temp(:,1:floor(size(v_temp,2)/(tr*1000)));
                                    %      v_tr_temp=resample(v_temp',1,tr*1000)';
                                    %     v_tr_temp=v_tr_temp(:,1:floor(size(v_temp,2)/(tr*1000)));
                                    
                                    if simi==1;
                                        tn=size(v_tr_temp,2)-1;
                                        v_hrf_tr=v_hrf_tr_temp(:,1:tn);
                                        v_tr=v_tr_temp(:,1:tn);
                                        v=v_temp(:,1:(tn*tr*1000));
                                        v_hrf=v_hrf_temp(:,1:(tn*tr*1000));
                                        onsetsVectors=onsetsVectors_temp(:,1:(tn*tr*1000));
                                    else
                                        tn=min(size(v_tr_temp,2)-1,size(v_tr,2));
                                        v_hrf_tr(:,1:tn,simi)=v_hrf_tr_temp(:,1:tn);
                                        v_tr(:,1:(tn),simi)=v_tr_temp(:,1:(tn));
                                        v(:,1:(tn*tr*1000),simi)=v_temp(:,1:(tn*tr*1000));
                                        v_hrf(:,1:(tn*tr*1000),simi)=v_hrf_temp(:,1:(tn*tr*1000));
                                        onsetsVectors(:,1:(tn*tr*1000),simi)=onsetsVectors_temp(:,1:(tn*tr*1000));
                                    end
                                end
                                if simN==1;
                                    save(sprintf('%s/simulation/sim%03d/pause%.2f_pauseL%.1f_m%d_v%.2f_sr%.2f_%s.mat',expdir,simN,pauseEffect,pauseLen,u,var,speechR, accumulationF),'v_hrf_tr','v_tr','v','v_hrf');
                                else
                                    save(sprintf('%s/simulation/sim%03d/pause%.2f_pauseL%.1f_m%d_v%.2f_sr%.2f_%s.mat',expdir,simN,pauseEffect,pauseLen,u,var,speechR, accumulationF),'v_hrf_tr','v_tr');
                                end
                        
                        end
                    end
                end
            end
        end
    end
end


