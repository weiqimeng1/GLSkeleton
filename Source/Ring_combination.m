function [Combined_ring, flag]=Ring_combination(main_ring, collapsed_ring, mainPInd, deletedPInd, pts, normalvec)
    if mainPInd==deletedPInd
        warning('no point to be deleted')
%         disp('current main point ind & main ring')
%         disp(mainPInd)
%         disp(main_ring)
%         disp('current main point ind & current collapsed ring')
%         disp(deletedPInd)
%         disp(collapsed_ring)
%         flag=0;
%         Combined_ring=main_ring;
        return;
    end
    Unioned_ring=union(main_ring, collapsed_ring);
    %*********** Point Candidate ************************************%
    temp_ind=setdiff(Unioned_ring,deletedPInd);
    temp_ind=[mainPInd,setdiff(temp_ind,mainPInd)];
    %************ Candidate coordinates*******************************%
    if length(temp_ind)<3
        flag=0;
        Combined_ring=main_ring;
        return
    end
    temp_points=pts(temp_ind,:);


    %************ inspiration from compute_point_point_ring***********%
    m=norm(normalvec);
    if m<eps
        coefs = pca(temp_points);
        x = [temp_points * coefs(:, 1), temp_points * coefs(:, 2)];
    else
        normalvec=normalvec';
        if normalvec(3)<eps %eps are numerical errors in Matlab, consider it as 0
            basevec1=[normalvec(2);-normalvec(1);0]/norm([normalvec(2),-normalvec(1),0]);
            vec2=cross(normalvec,basevec1);
            basevec2=vec2/norm(vec2);
        else
            vec1=[1;1;-(normalvec(1)+normalvec(2))/normalvec(3)];
            vec2=cross(normalvec,vec1);
            basevec1=vec1/norm(vec1);basevec2=vec2/norm(vec2);
        end
        x = [temp_points * basevec1, temp_points * basevec2];
    end
    TRI = delaunayn(x, {'QJ','Pp'});
    [row,~] = find(TRI == 1);
    % triangles connected centre
    temp = TRI(row,:);
    
    temp = sort(temp,2);
    temp = temp(:,2:end);
    
    x=temp(:);
    x=sort(x);

    d=diff([x;max(x)+1]);
    count = diff(find([1;d]));
    y =[x(d&1) count];
    n_sorted_index = size(y,1);
    start = find(count==1);
    if ~isempty(start) 
        want_to_find = y(start(1),1);
    else
        want_to_find = temp(1,1);
        n_sorted_index = n_sorted_index+1; 
    end
    
    j = 0;    
    sorted_index = zeros(1,n_sorted_index);
    while j < n_sorted_index
        j = j+1;
        sorted_index(j) = want_to_find;
        [row,col] = find(temp == want_to_find);
        if ~isempty(col)
            if col(1) == 1
                want_to_find = temp(row(1),2);
                temp(row(1),2) = -1;
            else
                want_to_find = temp(row(1),1);
                temp(row(1),1) = -1;
            end    
        end
    end
    Combined_ring = temp_ind(sorted_index);
    if length(Combined_ring)<=3
        Combined_ring=main_ring;
        flag=0;
%         disp('ring points number less than 3')
        return
    else
        flag=1;
    end
end