function dist=Eucdist(C,V)
if size(C,1)==1
    vlen=size(V,1);
    temp=repmat(C,vlen,1)-V;
    temp1=sum(temp.*temp,2);
    temp2=sqrt(temp1);
    dist=temp2;
else
    temp=C-V;
    temp1=sum(temp.*temp,2);
    temp2=sqrt(temp1);
    dist=temp2;
end