%% Settings
close all;  clear; clc;
path('io',path);
path('KDtree',path);
path('Visualisation',path);
path('Source',path);
options.SHOW_CONTRACTION_PROGRESS=false;
options.USING_POINT_RING = GS.USING_POINT_RING;
options.iterate_time = 10;
%%
% Load point cloud
contraction_cost=zeros(2,1);
% PD=pcread("data\horse_v1987.ply");
PD=pcread('\path\to\point\cloud');

P.pts=double(PD.Location);
P.faces=[];

%% Noise
%********************************************add noise**********************
% Define our signal - a column vector.
a=P.pts;
% Make the spread of the Gaussians be 20% of the a values
sigmas = 0 * a; % Also a column vector.
% Create the noise values that we'll add to a.
[m,n]=size(a);
randomNoise = randn(m,n) .* sigmas;
% Add noise to a to make an output column vector.
P.pts = a + randomNoise;

%% Laplacian-based contraction with point cloud reduction
P.npts = size(P.pts,1);
P.radis = ones(P.npts,1);
P.pts = GS.normalize(P.pts);
[P.bbox, P.diameter, P.area] = GS.compute_bbox(P.pts);
P.k_knn = GS.compute_k_knn(P.npts);
P.rings = compute_point_point_ring(P.pts, P.k_knn, []);

[P, t, initWL, WC, sl, Point_label, reducNum] = laplacian_contraction_with_reduction(P, options);

figure(1),
movegui('northeast');axis off;axis equal;set(gcf,'Renderer','OpenGL'); view3d rot;hold on;set(gcf,'color','white');
camorbit(0,0,'camera'); axis vis3d; view(-90,0);
h1 = scatter3(P.pts(:,1),P.pts(:,2), P.pts(:,3),2,[0 0.2235 0.3705],'filled');
h2 = scatter3(P.cpts(:,1),P.cpts(:,2), P.cpts(:,3),10,[0.9290 0.6940 0.1250],'filled');
legend('Original', 'contracted')
hold off

%% Extracting curve skeleton from contracted points
P.sample_radius = P.diameter*0.02; % fathest point sampling radius
P = extract_curve_from_contracted_points(P,P.sample_radius, 1); 

%% Visaulisation
showoptions.colorp=[0.8500 0.3250 0.0980];showoptions.colore=[0 0.2235 0.3705];
showoptions.sizep=400;showoptions.sizee=4;
figure,
plot_skeleton(P.spls, P.spls_adj, showoptions);
axis off;axis equal;set(gcf,'Renderer','OpenGL');view(0,90);view3d rot;