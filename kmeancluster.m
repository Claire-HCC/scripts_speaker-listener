eva = evalclusters(r,'kmeans','silhouett','KList',[1:10],'Distance','correlation')
[idx,C,sumd,D] = kmeans(r,2,'Distance','correlation');


r_temp=[];
for ei=1:4;
    exp=experiments{ei};
    load(['Y:\claire\speaker-listener\' exp '\fmri\pattern_lagcorr\tr\roi\mor\SLg\SL_lag-10-10.mat'])
   r_temp=[r_temp; r];
    

end

    eva = evalclusters(r_temp,'kmeans','silhouett','KList',[1:10],'Distance','correlation')
    [idx,C,sumd,D] = kmeans(r_temp,4,'Distance','correlation');
    