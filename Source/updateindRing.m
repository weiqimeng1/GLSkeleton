function updated_ind_ring=updateindRing(ind_ring, Poind_ind)
    prior_npts=length(Poind_ind);
    tempA=zeros(1,prior_npts);
    ind=1;
    for i=1:prior_npts
        if ind_ring(i)>0
            tempA(i)=ind;
            ind=ind+1;
        else
            tempA(i)=ind_ring(i);
        end
    end
    for i=1:prior_npts
        if tempA(i)<0
%             disp(tempA(-tempA(i)));
            tempA(i)=tempA(-tempA(i));
        end
    end
    updated_ind_ring=tempA;
end