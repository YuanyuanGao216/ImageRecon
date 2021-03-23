clear all
close all
clc
%% add snirf path
Snirf_path = 'C:\Users\sy100\Dropbox (BOSTON UNIVERSITY)\My PC (DESKTOP-J88OQCN)\Documents\GitHub\snirf_homer3';
Old_folder = cd(Snirf_path);
setpaths
cd(Old_folder)

%% load snirf data
acquired_path = 'G:\My Drive\BOAS\Joe_data\finger_tapping_nirx_hd\David\';
file = '2020-12-25_006.snirf';

acquired = SnirfClass([acquired_path,file]);
data = acquired.data;
probe = acquired.probe;

mlActMan =  cell(1,1);
mlActMan{1,1} = ones(100,1);

tIncMan =  cell(1,1);
tIncMan{1,1} = ones(size(acquired.data.time));

stim = acquired.stim;
Aaux = [];
tIncAuto = [];
rcMap = [];

%% homer3 GLM processing
Homer3_path = 'C:\Users\sy100\Dropbox (BOSTON UNIVERSITY)\My PC (DESKTOP-J88OQCN)\Documents\GitHub\Homer3';
Old_folder = cd(Homer3_path);
setpaths
cd(Old_folder)

mlActAuto = hmrR_PruneChannels(data,probe,mlActMan,tIncMan,[0  10000000],2,[0  45]);
pause(10)
dod = hmrR_Intensity2OD(data);
pause(0.1)
dod = hmrR_BandpassFilt(dod,0,0.5);
dc = hmrR_OD2Conc(dod,probe,[6  6  1]);
[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,dA,tbasis,tHRF] = hmrR_GLM_time_basis(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,[-2  20],1,3,[0.1 3.0 10.0 1.8 3.0 10.0],12,2,3);
% tbasis=zeros(ntHRF,nB,nConc);
%% Calculate L matrix
[Adot ,Adot_scalp, E] = Get_A_matrix(acquired_path);

%% Calculate T matrix
T_HbO_brain = squeeze(dA(:,:,1)); %4146*2
T_HbR_brain = squeeze(dA(:,:,2)); %4146*2

fs = 1/(dc.time(2)-dc.time(1));
frq1=(1:19)/20/fs;
T_scalp = zeros(length(dc.time),2*length(frq1));

for idx=1:length(frq1)
    T_scalp(:,1+(idx-1)*2)=sin(dc.time*2*pi*frq1(idx))';
    T_scalp(:,2+(idx-1)*2)=cos(dc.time*2*pi*frq1(idx))';
end

%% Calculate G matrix
At_file = 'atlasViewer.mat';
atlasViewer = load([acquired_path,At_file]);

brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
scalp_vertices = atlasViewer.fwmodel.mesh_scalp.vertices; %9563*3

% determine kernel positions
threshod = 30;
threshod_scalp = 120;
[brain_vertices_new] = down_sample_vertices(brain_vertices, threshod);
[scalp_vertices_new] = down_sample_vertices(scalp_vertices, threshod_scalp);

% visualize the new vertices on old mesh, trimesh
figure('name','brain')
T = atlasViewer.fwmodel.mesh.faces;
trimesh(T,brain_vertices(:,1),brain_vertices(:,2),brain_vertices(:,3))
hold on
plot3(brain_vertices_new(:,1),brain_vertices_new(:,2),brain_vertices_new(:,3),'ro');

figure('name','scalp')
T = atlasViewer.fwmodel.mesh_scalp.faces;
trimesh(T,scalp_vertices(:,1),scalp_vertices(:,2),scalp_vertices(:,3))
hold on
plot3(scalp_vertices_new(:,1),scalp_vertices_new(:,2),scalp_vertices_new(:,3),'ro');

% make kernels
sigma = 30;
sigma_scalp = 120;
G = make_kernel_matrix(brain_vertices_new, brain_vertices,sigma);
G_scalp = make_kernel_matrix(scalp_vertices_new, scalp_vertices,sigma_scalp);

%% H matrix
% Now we have:
% T_HbO_brain 4146,2
% T_HbR_brain 4146,2
% T_scalp 4146, 38
% G 122,20004
% G_scalp 13, 9562
% Adot 50,20004,2, 
% Adot_scalp 50, 9562, 2, 
% E 2,2
% OD 4146, 100
OD = dod.dataTimeSeries;
n_channels = size(OD,2);
time_points = size(OD,1);
Amatrix_brain = [squeeze(Adot(:,:,1))*E(1,1) squeeze(Adot(:,:,1))*E(1,2);
            squeeze(Adot(:,:,2))*E(2,1) squeeze(Adot(:,:,2))*E(2,2)]; % 100, 59132
Amatrix_scalp = [squeeze(Adot_scalp(:,:,1))*E(1,1) squeeze(Adot_scalp(:,:,1))*E(1,2);
            squeeze(Adot_scalp(:,:,2))*E(2,1) squeeze(Adot_scalp(:,:,2))*E(2,2)]; % 100, 59132

n_g = size(G,1);
n_t = size(T_HbO_brain,2);
% H_brain = zeros(n_channels*time_points,n_g*n_t);
% j = 1;
% for i_g = 1:n_g
%     fprintf('%d/%d\n',i_g,n_g)
%     for i_t = 1:n_t
%         HbO = kron(G(i_g,:),T_HbO_brain(:,i_t));
%         HbR = kron(G(i_g,:),T_HbR_brain(:,i_t));
%         gt = [HbO,HbR]';
%         Agt = Amatrix_brain*gt;
%         H_brain(:,j) = reshape(Agt,[],1);
%         j = j + 1;
%     end
% end
% 
% n_g = size(G_scalp,1);
% n_t = size(T_scalp,2);
% H_scalp = zeros(n_channels*time_points,n_g*n_t);
% j = 1;
% for i_g = 1:n_g
%     fprintf('%d/%d\n',i_g,n_g)
%     for i_t = 1:n_t
%         HbO = kron(G_scalp(i_g,:),T_scalp(:,i_t));
%         HbR = kron(G_scalp(i_g,:),T_scalp(:,i_t));
%         gt = [HbO,HbR]';
%         Agt = Amatrix_scalp*gt;
%         H_scalp(:,j) = reshape(Agt,[],1);
%         j = j + 1;
%     end
% end
%% H matrix is too large, we calculated H'*H line by line here

Y = reshape(OD',[],1);

HTH = zeros(size(G,1)*size(T_HbO_brain,2)+ size(G_scalp,1)*size(T_scalp,2));
HTY = 
j1 = 0;

for i_g = 1:(size(G,1)+ size(G_scalp,1))
    if i_g <= size(G,1)
        g = G(i_g,:);
        n_t = size(T_HbO_brain,2);
        T_HbO = T_HbO_brain;
        T_HbR = T_HbR_brain;
        Amatrix = Amatrix_brain;
    else
        g = G_scalp(i_g - size(G,1),:);
        n_t = size(T_scalp,2);
        T_HbO = T_scalp;
        T_HbR = T_scalp;
        Amatrix = Amatrix_scalp;
    end
    
    for i_t = 1:n_t
        HbO = kron(g,T_HbO(:,i_t));
        HbR = kron(g,T_HbR(:,i_t));
        gt = [HbO,HbR]';
        Agt = Amatrix*gt;
        H1 = reshape(Agt,[],1);
        j1 = j1+1;
        HTY = H1'*Y;
        % loop all cols of H
        j2 = 0;
        for i_g2 = 1:(size(G,1)+ size(G_scalp,1))
            if i_g2 <= size(G,1)
                g2 = G(i_g,:);
                n_t2 = size(T_HbO_brain,2);
                T_HbO2 = T_HbO_brain;
                T_HbR2 = T_HbR_brain;
                Amatrix2 = Amatrix_brain;
            else
                g2 = G_scalp(i_g - size(G,1),:);
                n_t2 = size(T_scalp,2);
                T_HbO2 = T_scalp;
                T_HbR2 = T_scalp;
                Amatrix2 = Amatrix_scalp;
            end
            for i_t2 = 1:n_t2
                HbO2 = kron(g,T_HbO2(:,i_t));
                HbR2 = kron(g,T_HbR2(:,i_t));
                gt2 = [HbO2,HbR2]';
                Agt2 = Amatrix2*gt;
                H2 = reshape(Agt2,[],1);
                j2 = j2+1;
                
                HTH(j1, j2) = H1'*H2;
                HTY
            end
        end
    end
end
%% do inverse
H = [H_brain,H_scalp];



beta = (H'*Y)'/(H'*H); 
% regularization [lamda_brain*I 0;0 lamda_scalp*I]
%% approximate the conc
n_g = size(G,1);
n_t = size(T_HbO_brain,2);
beta_brain = beta(1:n_g*n_t);

Conc_matrix = zeros(size(G,2)*2,size(tbasis,1));

j = 1;

for i_g = 1:n_g
    fprintf('%d/%d\n',i_g,n_g)
    for i_t = 1:n_t
        HbO = kron(G(i_g,:),tbasis(:,i_t,1));
        HbR = kron(G(i_g,:),tbasis(:,i_t,2));
        gt = [HbO,HbR]';
        Conc_matrix = Conc_matrix + gt.*beta_brain(j);
        j = j + 1;
    end
end

intensity_HbO = Conc_matrix(1:size(G,2),:);
intensity_HbR = Conc_matrix(size(G,2)+1:end,:);

save('intensity.mat','intensity_HbO','intensity_HbR');
%% visualization
Homer3_path = 'C:\Users\sy100\Dropbox (BOSTON UNIVERSITY)\My PC (DESKTOP-J88OQCN)\Documents\GitHub\Homer3';
Old_folder = cd(Homer3_path);
setpaths(0)
cd(Old_folder)

AtlasViewer_path = 'C:\Users\sy100\Dropbox (BOSTON UNIVERSITY)\My PC (DESKTOP-J88OQCN)\Documents\GitHub\AtlasViewer';
Old_folder = cd(AtlasViewer_path);
setpaths
cd(Old_folder)

mesh = atlasViewer.imgrecon.mesh;
light_onoff = 'on';
visible = 'on';
axes_order = [1 2 3];


close all
figure
subplot(321)
[h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>-1,1)), light_onoff, visible, axes_order);
colorbar
% caxis([-8 100])
title('-1s')

subplot(322)
[h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>0,1)), light_onoff, visible, axes_order);
colorbar
% caxis([-8 100])
title('0s')

subplot(323)
[h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>5,1)), light_onoff, visible, axes_order);
colorbar
% caxis([-8 100])
title('5s')

subplot(324)
[h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>10,1)), light_onoff, visible, axes_order);
colorbar
% caxis([-8 100])
title('10s')

subplot(325)
[h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>15,1)), light_onoff, visible, axes_order);
colorbar
% caxis([-8 100])
title('15s')

subplot(326)
[h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>20,1)), light_onoff, visible, axes_order);
colorbar
% caxis([-8 100])
title('20s')

set(gcf,'position',[ 1001         135        1147        1204])