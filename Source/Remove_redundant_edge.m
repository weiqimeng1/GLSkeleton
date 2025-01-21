function [new_spls,new_adj]=Remove_redundant_edge(spls,adj)
npts=size(spls,1);
new_adj=adj;
new_spls=spls;
new_adj=new_adj-diag(ones(npts,1));
for i=1:npts
    nlist=new_adj(i,:);
    nlist=nlist(nlist>0);
    if length(nlist)==2
        v1=spls(i,:)-spls(nlist(1),:);
        v2=spls(i,:)-spls(nlist(2),:);
        cosAngle=dot(v1,v2)/(norm(v1)*norm(v2));
        if cosAngle<-0.98
            new_spls(i,:)=nan;
            new_adj(i,nlist(1))=0;
            new_adj(i,nlist(2))=0;
            new_adj(nlist(1),i)=0;
            new_adj(nlist(2),i)=0;
            new_adj(nlist(2),nlist(1))=1;
            new_adj(nlist(1),nlist(2))=1;
        end
    end
end
new_adj=new_adj+diag(ones(npts,1));
end