function [newP, newCpts, flag, Point_ind] = Remove_points_in_rings(P, cpts, Checklist, MDI,MDV, GammaVec)
    newP=P;
    newCpts=cpts;
    Point_ind=ones(newP.npts,1);
    Combinedone_list=zeros(newP.npts,1);
    Ind_ring=1:P.npts;
    %**********************if there are rings have been collapsed*************%
    collapsed_flag=0;
    %********************Get collapsed Point index************************
    for i=1:newP.npts
        if Combinedone_list(i)==1
            continue;
        end
        cuvdiffVec=GammaVec(1,:)-sum(GammaVec(newP.rings{i},:),1)./length(newP.rings{i});
        check_value=Checklist(i)+3*Vec_Square(cuvdiffVec);
        CombineI=newP.rings{i}(MDI(i));
        if check_value<1e5*(MDV(i))^2 && Combinedone_list(CombineI)~=1
            [new_ring, crflag]=Ring_combination(newP.rings{i},newP.rings{CombineI}, i, CombineI, newP.pts,GammaVec(i,:));
            if crflag
                newP.rings{i}=new_ring;
                collapsed_flag=1;
                Point_ind(CombineI)=0;
                Ind_ring(CombineI)=-Ind_ring(i);
                Combinedone_list(i)=1;
                Combinedone_list(CombineI)=1;
            end
        end
    end
    %**********************update points& rings***********************%
    if collapsed_flag
        newP.pts=newP.pts(Point_ind&1,:); newCpts=newCpts(Point_ind&1,:);

        newP.npts=size(newP.pts,1);
        %************urflag: if the rings have been updated correctly*******************%
       
        [newP.rings, urflag]=update_rings_after_manipulation(newP.rings, Point_ind, Ind_ring);
        if ~urflag
            warning('urflag wrong')
            newP=P;newCpts=cpts;Point_ind=1:P.npts; flag=false;return;
        end
    else
        % newP=P is not needed
        flag=false;
        warning('no collapse')
        return;
    end
    flag=true;
end

function Normlist=Vec_Square(Veclist)
Normlist=sum(Veclist.*Veclist,2);
end