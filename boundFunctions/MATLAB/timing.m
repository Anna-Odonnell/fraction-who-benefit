rng(3572032359)

nLevels = 20;
sampSize = 5000;

YT = randsample(nLevels,sampSize,true);

YC = randsample(nLevels,sampSize,true);

tic
[~,~,l,u,eps] = boundsNoCov_res(YT, YC, 10^5, 10^5);
toc