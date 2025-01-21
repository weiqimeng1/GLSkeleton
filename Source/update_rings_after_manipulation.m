function [newrings, flag]=update_rings_after_manipulation(rings, Poind_ind, ind_ring)
    %***************Check if update is required****************%
    if all(Poind_ind)
        flag=0;
        return
    end

    %**************update rings*****************************%
    flag=1;
    newrings=rings(Poind_ind&1);
    updated_ind_ring=updateindRing(ind_ring, Poind_ind);
    for i=1:sum(Poind_ind)
        newrings{i}=updated_ind_ring(newrings{i});
        ring_differ=diff([newrings{i},newrings{i}(end)+1]);
        newrings{i}=newrings{i}(logical(ring_differ));
    end
end