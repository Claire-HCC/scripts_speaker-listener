
close all
clear all
% loc='cluster';
set_parameters;

exp='sherlock';
tb=readtable([expdir exp '/sound/transcripts_'  exp '_Claire.xlsx']);

% get hrf
% xBF.dt=tr(ei);
% xBF.name='hrf (with time and dispersion derivatives)';
% bf = spm_get_bf(xBF);
% tr=1.5
load([expdir '/simulation/bf.mat']);

tr=1.5;
simN=1;
NRFs={'linear', 'log', 'square','time'};
pauses_sec=[3 1.5 0];
wdN=3000;
scramble=0;

for psi=1%:length(pauses_sec);
    pause_sec=pauses_sec(psi);
    
    for u=3;%2:4
        for vi=1.4;%[1.6 1.9];%[0.1:0.3:1.5];
            for speechR=1;%[0.5 1 2]; % 1 is the speech rate of sherlock;
                for af=4%:3;
                    NRF=NRFs{af};
                    if exist(sprintf('%s/simulation/sim%03d/pause%.1f_m%d_v%.2f_sr%.2f_%s.mat',expdir,simN,pause_sec,u,vi,speechR, NRF))~=0;
                        load(sprintf('%s/simulation/sim%03d/pause%.1f_m%d_v%.2f_sr%.2f_%s.mat',expdir,simN,pause_sec,u,vi,speechR, NRF),'onsetsVectors','v','v_hrf_sec');
                        
                        % remove pause
                        %    onsetsVectors(:,v_temp(levn,:)==0)=[];
                        
                        for simi=1:simN;
                            onsetsVectors_temp=onsetsVectors(:,:,simi);
                            v_temp=v(:,:,simi);
                            
                            [levn,tn]=size(onsetsVectors_temp);
                            
                            onsetsVectors_scrambled_temp=zeros(size(onsetsVectors_temp));
                            v_scrambled_temp=zeros(size(v_temp));
                            
                            for lev=1:(levn);
                                
                                onsets=[1 find(onsetsVectors_temp(lev,:)==1)];
                                offsets=[find(onsetsVectors_temp(lev,:)==1)-1 tn];
                                n=length(onsets);
                                scrambledi=randperm(n);
                                onsets_scrambled=onsets(scrambledi);
                                offsets_scrambled=offsets(scrambledi);
                                
                                otemp=[];
                                vtemp=[];
                                for ni=1:n;
                                    idx=onsets_scrambled(ni):offsets_scrambled(ni);
                                    otemp=[otemp onsetsVectors_temp(lev,idx)];
                                    vtemp=[vtemp v_temp(lev,idx)];
                                end
                                
                                onsetsVectors_scrambled_temp(lev,:)=otemp;
                                v_scrambled_temp(lev,:)=vtemp;
                            end
                            
                            % hrf convolution
                            v_hrf_scrambled=[];
                            for xi=1:size(v_scrambled_temp,1);
                                temp=conv(v_scrambled_temp(xi,:),bf.bf(:,1));
                                v_hrf_scrambled(xi,:)=temp;
                            end
                            
                            % downsampling to  1 sec resolution
                            v_hrf_sec_scrambled_temp=resample(v_hrf_scrambled',1,tr*1000)';
                            
                            v_hrf_sec_scrambled_temp=v_hrf_sec_scrambled_temp(:,1:floor(size(v_scrambled_temp,2)/(tr*1000)));
                            
                            if simi==1;
                                tn=size(v_hrf_sec_scrambled_temp,2);
                                v_hrf_sec_scrambled=v_hrf_sec_scrambled_temp;
                                v_scrambled=v_scrambled_temp(:,1:(tn*tr*1000));
                                onsetsVectors_scrambled=onsetsVectors_scrambled_temp(:,1:(tn*tr*1000));
                            else
                                tn=min(size(v_hrf_sec_scrambled_temp,2),size(v_hrf_sec_scrambled,2));
                                v_hrf_sec_scrambled(:,1:tn,simi)=v_hrf_sec_scrambled_temp(:,1:tn);
                                v_scrambled(:,1:(tn*tr*1000),simi)=v_scrambled_temp(:,1:(tn*tr*1000));
                                onsetsVectors_scrambled(:,1:(tn*tr*1000),simi)=onsetsVectors_scrambled_temp(:,1:(tn*tr*1000));
                            end
                        end
                        v_hrf_sec=v_hrf_sec_scrambled;
                        v=v_scrambled;
                        onsetsVectors=onsetsVectors_scrambled;
                        
                        save(sprintf('%s/simulation/sim%03d_scrambled/pause%.1f_m%d_v%.2f_sr%.2f_%s.mat',expdir,simN,pause_sec,u,vi,speechR, NRF),'v_hrf_sec','v','onsetsVectors');
                    end
                end
            end
        end
    end
end



