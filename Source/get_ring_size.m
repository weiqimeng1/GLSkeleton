function [ms, mi] = get_ring_size(pts,rings)
n = size(pts,1);
ms = zeros(n,1);
mi = zeros(n,1);
parfor i = 1:n
    ring = rings{i};    
    tmp = repmat(pts(i,:), length(ring),1) - pts(ring,:);
    [ms(i),mi(i)] = min( sum(tmp.^2,2).^0.5);
end