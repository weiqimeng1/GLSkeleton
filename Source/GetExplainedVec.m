function Evec=GetExplainedVec(PtCloud)
% PtCloud n*3
M=cov(PtCloud);
[~, E]=eig(M);
Evec=diag(E);