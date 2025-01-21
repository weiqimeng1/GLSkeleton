function [P, t, initWL, WC, sl, Point_label, reducNum] = laplacian_contraction_with_reduction(P,options,SHOW_CONTRACTION_PROGRESS)
% Laplacian-based contraction with point cloud reduction
% refer to Point Cloud Skeletons via Laplacian-Based Contraction
%
% inputs:
%   P.pts
%   P.faces
%   P.npts
%   P.k_knn: k of knn
%   P.rings:
%   options.WL = 3;% initial contraction weight
%   options.WC = 1;% initial attraction weight
%   options.sl: scalar for increasing WL in each iteration
%   options.tc: Termination Conditions for total area ratio
%   options.iterate_time = 10; max iteration steps
%
% outputs:
%   cpts: contracted vertices
%   Point_label: Point removed in iterations and remained are marked with
%   different number
%   reducNum: number of removed points
%
% notes:
%   Reference: https://github.com/taiya/cloudcontr
%   Adapted by: Q. Wen
%   Wen, Qingmeng, et al. "GLSkeleton: A geometric Laplacian-based skeletonisation framework for object point clouds." IEEE Robotics and Automation Letters (2024).


if nargin < 1
    clear;clc;close all;
    P.filename = '../data/simplejoint_v4770.off';
    options.USING_POINT_RING = GS.USING_POINT_RING;
    options.iterate_time = 10;
    [P.pts,P.faces] = read_mesh(P.filename);
    P.npts = size(P.pts,1);
    P.pts = GS.normalize(P.pts);
    [P.bbox, P.diameter] = GS.compute_bbox(P.pts);

    P.k_knn = GS.compute_k_knn(P.npts);
    atria = nn_prepare(P.pts);
    [P.knn_idx, P.knn_dist] = nn_search(P.pts, atria, P.pts, P.k_knn);
    P.rings = compute_point_point_ring(P.pts,P.k_knn, P.knn_idx);
end

%##########################################################################
%% settings
%##########################################################################
orginal_npts=P.npts;
ctr=1;
% visual debug conditions
RING_SIZE_TYPE = 1;%1:min, 2:mean, 3:max
Laplace_type = 'conformal';%conformal%combinatorial%spring%mvc
cflag=true;

% setting
tc = getoptions(options, 'tc', GS.CONTRACT_TERMINATION_CONDITION);
iterate_time = getoptions(options, 'iterate_time', GS.MAX_CONTRACT_NUM);

initWL = getoptions(options, 'WL', GS.compute_init_laplacian_constraint_weight(P,Laplace_type));
% set initial attraction weights according to different type of discrete
% Laplacian
WC = getoptions(options, 'WC', 1);

WH = ones(P.npts, 1)*WC; % initialise weights
sl = getoptions(options, 'sl', GS.LAPLACIAN_CONSTRAINT_SCALE); % scale factor for WL in each iteration! in original paper is 2;
WL = initWL;%*sl;

%% init iteration
t = 1; % current iteration step
%% left side of the equation
L = -compute_point_laplacian(P.pts,Laplace_type,P.rings, options);%conformal;%spring
A = [L.*WL;sparse(1:P.npts,1:P.npts, WH)];
% right side of the equation
b = [zeros(P.npts,3);sparse(1:P.npts,1:P.npts, WH)*P.pts];
cpts = (A'*A)\(A'*b);
EVec=GetExplainedVec(cpts);
e2=(2*EVec(3)+EVec(2)+EVec(1))/(sum(EVec)*2);
Point_label=ones(P.npts,1);
ind_remained=1:P.npts;
reducNum=0;

if options.SHOW_CONTRACTION_PROGRESS
    figure;movegui('northeast');axis off;axis equal;set(gcf,'Renderer','OpenGL'); view3d rot;hold on;set(gcf,'color','white');
    camorbit(0,0,'camera'); axis vis3d; view(-90,0);
    h1 = scatter3(P.pts(:,1),P.pts(:,2), P.pts(:,3),10,'b','filled');
    h2 = scatter3(cpts(:,1),cpts(:,2), cpts(:,3),10,'r','filled');
    title(['iterate ',num2str(t),' time(s)'])
end
%%
sizes = GS.one_ring_size(P.pts, P.rings, RING_SIZE_TYPE);   % min radius of 1-ring
[size_new, MDI] = get_ring_size(cpts, P.rings);
a(t) = sum(size_new)/sum(sizes);

%%%%%%%%%%%%%%%%%%%%%%Point reduction%%%%%%%%%%%%%%%%%%%%%%%%%%
if cflag
    disvec=P.pts-cpts;
    CGamma_vec=L*P.pts;
    Controled_ratio=1.2*a(end);
    S=size_new-Controled_ratio*sizes;
    S(S<0)=0.0001;
    SizeRatio=S./sizes;
    Checklist=Vec_Square(sl^2*initWL.*(CGamma_vec.*SizeRatio)) ...
        -Vec_Square(1./SizeRatio.*((Controled_ratio*sizes./(sizes-size_new).*disvec)));
    [P, cpts, cflag, pind] = Remove_points_in_rings(P, cpts, Checklist,MDI,size_new, CGamma_vec);
    reducNum=reducNum+1;
    Point_label(ind_remained(pind&1))=reducNum;
    ind_remained=ind_remained(pind&1);
    ctc=P.npts/orginal_npts;
    if cflag
        size_new=size_new(pind&1);
        sizes=sizes(pind&1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while t<iterate_time
    L = -compute_point_laplacian(cpts,Laplace_type, P.rings, options);%conformal

    WL = sl*WL;
    if WL>GS.MAX_LAPLACIAN_CONSTRAINT_WEIGHT
        WL=GS.MAX_LAPLACIAN_CONSTRAINT_WEIGHT;
    end % from Oscar08's implementation, 2048
    WH = WC.*(sizes./size_new); % update attraction weights

    WH(WH>GS.MAX_POSITION_CONSTRAINT_WEIGHT) = GS.MAX_POSITION_CONSTRAINT_WEIGHT;% from Oscar08's implementation, 10000

    A = real([WL*L;sparse(1:P.npts,1:P.npts, WH)]);

    % update right side of the equation
    b = [zeros(P.npts,3);sparse(1:P.npts,1:P.npts, WH)*cpts];
    tmp = (A'*A)\(A'*b);

    [size_new, MDI] = get_ring_size(tmp, P.rings);
    a(end+1) = sum(size_new)/sum(sizes);

    tmpbox = GS.compute_bbox(tmp);
    if sum( (tmpbox(4:6)-tmpbox(1:3))> ((P.bbox(4:6)-P.bbox(1:3))*1.2) ) > 0
        break;
    end

    if isnan(a(end))
        break;
    else
        %%%%%%%%%%%%%%%%break condition%%%%%%%%%%%%%%%%%%%%%%%%
        e1=e2;
        EVec=GetExplainedVec(cpts);
        e2=(2*EVec(3)+EVec(2)+EVec(1))/(sum(EVec)*2);
        explained_diff=abs((e2-e1)*(1-ctr+ctc));
        if t>=3 && explained_diff<0.02
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cpts_prior=cpts;
        cpts = tmp;
        %%%%%%%%%%%%%%%%%%%%Point reduction%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if cflag
            disvec=P.pts-cpts;
            CGamma_vec=L*cpts_prior;
            Controled_ratio=1.2*(a(end-1)-a(end));
            S=size_new-Controled_ratio*sizes;
            S(S<0)=0.0001;
            SizeRatio=S./sizes;
            Checklist=Vec_Square(sl^(t+2)*initWL.*(CGamma_vec.*SizeRatio)) ...
                -Vec_Square(1./SizeRatio.*((Controled_ratio*sizes./(sizes-size_new).*disvec)));
            [P, cpts, cflag, pind] = Remove_points_in_rings(P, cpts, Checklist, MDI,size_new, CGamma_vec);
            reducNum=reducNum+1;
            Point_label(ind_remained(pind&1))=reducNum;
            ind_remained=ind_remained(pind&1);
            ctr=ctc;
            ctc=P.npts/orginal_npts;
            if cflag
                size_new=size_new(pind&1);
                sizes=sizes(pind&1);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    t = t+1;

    if options.SHOW_CONTRACTION_PROGRESS
        % Show point cloud change

        delete(h1);delete(h2);
        h1 = scatter3(P.pts(:,1),P.pts(:,2), P.pts(:,3),10,'b','filled');
        h2 = scatter3(cpts(:,1),cpts(:,2), cpts(:,3),10,'r','filled');
        %legend('orignal points','contracted points');
        title(['iterate ',num2str(t),' time(s)']); drawnow;
    end
end
clear tmp;

if options.SHOW_CONTRACTION_PROGRESS
    figure;
    plot(a);xlabel('Iteration times');ylabel('Ratio of original and current volume');
end
P.cpts=cpts;
end


function Normlist=Vec_Square(Veclist)
Normlist=sum(Veclist.*Veclist,2);
end