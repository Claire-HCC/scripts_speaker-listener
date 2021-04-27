
us=3; % 2:4; 3000 words are often not enough to create level 6 with more than to 2 consituents.
vars=0.5;
speechRs=1;
accumulationFs={'linear'}; % {'linear', 'log','expn', 'triangle','time'};
pauseLens=3;
pauseEffects=0.1; % activity at pause period = min acitivity-sd of activity *pauseEffect
simulation_pallier(us,vars,[0.5 1 1.5],accumulationFs,pauseLens,pauseEffects);
simulation_pallier_circularlagcorr;
simulation_pallier_stats;
simulation_pallier_circularlagcorr_imagesc;
