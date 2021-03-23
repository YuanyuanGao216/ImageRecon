clear all
close all
clc

%% tool folder manage
add_snirf_path('C:\Users\sy100\Dropbox (BOSTON UNIVERSITY)\My PC (DESKTOP-J88OQCN)\Documents\GitHub\snirf_homer3')
add_homer_path()
% add_av_path()

%% add simulation fusion path
addpath(genpath('SIMULATION_FUSION'))
%% add snirf path
Snirf_path = 'C:\Users\sy100\Dropbox (BOSTON UNIVERSITY)\My PC (DESKTOP-J88OQCN)\Documents\GitHub\snirf_homer3';
Old_folder = cd(Snirf_path);
setpaths
cd(Old_folder)

%% load snirf data
acquired_path = ['data',filesep];
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
dod = hmrR_Intensity2OD(data);
dod = hmrR_BandpassFilt(dod,0,0.5);
dc = hmrR_OD2Conc(dod,probe,[6  6  1]);
[dcAvg,dcAvgStd,nTrials,dcNew,dcResid,dcSum2,beta,R,dA,tbasis,tHRF] = hmrR_GLM_time_basis(dc,stim,probe,mlActAuto,Aaux,tIncAuto,rcMap,[-2  20],1,3,[0.1 3.0 10.0 1.8 3.0 10.0],12,2,3);
% tbasis=zeros(ntHRF,nB,nConc);
%% Calculate L matrix
[Adot ,Adot_scalp, E] = Get_A_matrix(acquired_path);

%% Calculate T matrix
T_HbO_brain = squeeze(dA(:,1,1)); %4146*1 only the slow one
T_HbR_brain = squeeze(dA(:,1,2)); %4146*1

fs = 1/(dc.time(2)-dc.time(1));
frq1=(1:2)/2/fs;
T_scalp = zeros(length(dc.time),2*length(frq1));

for idx=1:length(frq1)
    T_scalp(:,1+(idx-1)*2)=sin(dc.time*2*pi*frq1(idx))';
    T_scalp(:,2+(idx-1)*2)=cos(dc.time*2*pi*frq1(idx))';
end


% time = data.time;
% stim_mark = stim(1,1).data(:,1);
% 
% stim_data = zeros(size(time));
% [sharedvals,idx] = intersect(time,stim_mark,'stable');
% stim_data(idx) = 1;
% 
% [To_HbO To_HbR] =  createtemporalbase(stim_data);
% [TO] = combineBasis(To_HbO,To_HbR);

%% Calculate G matrix
At_file = 'atlasViewer.mat';
atlasViewer = load([acquired_path,At_file]);

brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
scalp_vertices = atlasViewer.fwmodel.mesh_scalp.vertices; %9563*3

%% making mask here, output is mask and mask_scalp
figure
intensity = log10(sum(Adot(:,:,1),1));
faces = atlasViewer.fwmodel.mesh.faces;
trisurf(faces,brain_vertices(:,1),brain_vertices(:,2),brain_vertices(:,3),intensity,'edgecolor','none')
mask = find(intensity>-2);
fprintf('length is %d\n', length(mask))

figure
intensity_scalp = log10(sum(Adot_scalp(:,:,1),1));
faces_scalp = atlasViewer.fwmodel.mesh_scalp.faces;
trisurf(faces_scalp,scalp_vertices(:,1),scalp_vertices(:,2),scalp_vertices(:,3),intensity_scalp,'edgecolor','none')
mask_scalp = find(intensity_scalp>-2);
fprintf('scalp length is %d\n', length(mask_scalp))


% here i want to show what is left after mask
% figure
% trisurf(T,brain_vertices(:,1),brain_vertices(:,2),brain_vertices(:,3),'edgecolor',[0.1 0.1 0.1],'edgecolor','none')
% hold on
% trisurf(T(),brain_vertices(mask,1),brain_vertices(mask,2),brain_vertices(mask,3),'edgecolor',[0.1 0.1 0.1],'edgecolor','none')
%% Amatrix after mask
Adot = Adot(:, mask, :);
Adot_scalp = Adot_scalp(:, mask_scalp, :);

Amatrix_brain = [squeeze(Adot(:,:,1))*E(1,1) squeeze(Adot(:,:,1))*E(1,2);
            squeeze(Adot(:,:,2))*E(2,1) squeeze(Adot(:,:,2))*E(2,2)]; % 100, 40008
Amatrix_scalp = [squeeze(Adot_scalp(:,:,1))*E(1,1) squeeze(Adot_scalp(:,:,1))*E(1,2);
            squeeze(Adot_scalp(:,:,2))*E(2,1) squeeze(Adot_scalp(:,:,2))*E(2,2)]; % 100, 19124

%% G matrix after mask
% determine kernel positions
brain_vertices_masked = brain_vertices(mask,:);
scalp_vertices_masked = scalp_vertices(mask_scalp,:);

threshod = 15;
threshod_scalp = 80;
[brain_vertices_new] = down_sample_vertices(brain_vertices_masked, threshod);
[scalp_vertices_new] = down_sample_vertices(scalp_vertices_masked, threshod_scalp);

% visualize the new vertices on old mesh, trimesh
figure('name','brain')
faces = atlasViewer.fwmodel.mesh.faces;
trimesh(faces,brain_vertices(:,1),brain_vertices(:,2),brain_vertices(:,3))
hold on
plot3(brain_vertices_new(:,1),brain_vertices_new(:,2),brain_vertices_new(:,3),'ro');
set(gcf,'position',[10 10 560 420])

figure('name','scalp')
faces_scalp = atlasViewer.fwmodel.mesh_scalp.faces;
trimesh(faces_scalp, scalp_vertices(:,1),scalp_vertices(:,2),scalp_vertices(:,3))
hold on
plot3(scalp_vertices_new(:,1),scalp_vertices_new(:,2),scalp_vertices_new(:,3),'ro');
set(gcf,'position',[10 10 560 420])

% make kernels
sigma = 15;
sigma_scalp = 80;
G = make_kernel_matrix(brain_vertices_new, brain_vertices_masked,sigma);
G_scalp = make_kernel_matrix(scalp_vertices_new, scalp_vertices_masked,sigma_scalp);

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
%% H matrix is too large, we calculated H'*H line by line here

Y = reshape(OD',[],1); % Y = [414600, 1]

HTH = zeros(size(G,1)*size(T_HbO_brain,2)+ size(G_scalp,1)*size(T_scalp,2),'gpuArray'); % HTH = [240, 240]
HTY = zeros(size(G,1)*size(T_HbO_brain,2)+ size(G_scalp,1)*size(T_scalp,2), 1,'gpuArray'); % HTY = [240, 1]
j1 = 0;

tic
% ram was 1234
for i_g = 1:(size(G,1)+ size(G_scalp,1))
    fprintf('i_g is %d/%d\n', i_g, (size(G,1)+ size(G_scalp,1)))
    if i_g <= size(G,1)
        g = G(i_g,:); % G = [50, 1576], g = [1, 1576]
        n_t = size(T_HbO_brain,2);  
        T_HbO = T_HbO_brain; % T_HbO_brain = [4146,1],
        T_HbR = T_HbR_brain; % T_HbR_brain = [4146,1],
        Amatrix = Amatrix_brain; % Amatrix = [100, 3152]
    else
        g = G_scalp(i_g - size(G,1),:);% G_scalp = [19, 9562], g = [1, 9562]
        n_t = size(T_scalp,2);
        T_HbO = T_scalp;% T_HbO = [4146,38],
        T_HbR = T_scalp;% T_HbR = [4146,38],
        Amatrix = Amatrix_scalp;% Amatrix = [100, 19124]
    end
    
    g =  gpuArray(g);
    T_HbO = gpuArray(T_HbO);
    T_HbR = gpuArray(T_HbR);
    Amatrix = gpuArray(Amatrix);
    
    for i_t = 1:n_t
%         HbO = kron(g,T_HbO(:,i_t)); % HbO = [4146, 20004],
%         HbR = kron(g,T_HbR(:,i_t)); % HbR = [4146, 20004],
%         gt = [HbO,HbR]';
        gt = [kron(g,T_HbO(:,i_t)),kron(g,T_HbR(:,i_t))]'; % gt = [3152, 4146],
        Agt = Amatrix*gt; % Agt = [100, 4146]
        H1 = reshape(Agt,[],1);% H1 = [414600, 1]
        if isfile('H.mat')
            matObj.H(:,end+1) = gather(H1);
        else
            H = gather(H1);
            save('H.mat','H');
            matObj = matfile('H.mat','Writable',true);
        end
        
        j1 = j1+1;
        HTY(j1) = H1'*Y;
        % loop all cols of H
        j2 = 0;
        for i_g2 = 1:(size(G,1)+ size(G_scalp,1))
            if i_g2 <= size(G,1)
                g2 = G(i_g2,:); % g2 = [1, 1576]
                n_t2 = size(T_HbO_brain,2);
                T_HbO2 = T_HbO_brain; % T_HbO2 = [4146, 1]
                T_HbR2 = T_HbR_brain; % T_HbR2 = [4146, 1]
                Amatrix2 = Amatrix_brain; % Amatrix2 = [100, 3152]
            else
                g2 = G_scalp(i_g2 - size(G,1),:);
                n_t2 = size(T_scalp,2);
                T_HbO2 = T_scalp;
                T_HbR2 = T_scalp;
                Amatrix2 = Amatrix_scalp;
            end
            
            g2 =  gpuArray(g2);
            T_HbO2 = gpuArray(T_HbO2);
            T_HbR2 = gpuArray(T_HbR2);
            Amatrix2 = gpuArray(Amatrix2);
            
            for i_t2 = 1:n_t2
%                 HbO2 = kron(g2,T_HbO2(:,i_t2)); % HbO2 = [4146, 20004]
%                 HbR2 = kron(g2,T_HbR2(:,i_t2)); % HbR2 = [4146, 20004]
                gt2 = [kron(g2,T_HbO2(:,i_t2)),kron(g2,T_HbR2(:,i_t2))]'; % gt2 = [3152, 4146]
                Agt2 = Amatrix2*gt2; % Agt2 = [100, 4146]
                H2 = reshape(Agt2,[],1); % H2 = [414600, 1]
                j2 = j2+1;
                
                HTH(j1, j2) = H1'*H2;
            end
        end
    end
    toc
    
end

save('HTH.mat','HTH')
%% do inverse

beta = (HTY)'/(HTH); % beta [1, 240]
% regularization [lamda_brain*I 0;0 lamda_scalp*I]
%% approximate the conc
n_g = size(G,1);
n_t = size(T_HbO_brain,2);
beta_brain = beta(1:n_g*n_t);

Conc_matrix = zeros(size(G,2)*2,size(T_HbO_brain,1));

j = 1;

for i_g = 1:n_g
    for i_t = 1:n_t
        HbO = kron(G(i_g,:),T_HbO_brain(:,i_t));
        HbR = kron(G(i_g,:),T_HbR_brain(:,i_t));
        gt = [HbO,HbR]';
        Conc_matrix = Conc_matrix + gt.*beta_brain(j);
        j = j + 1;
    end
end

intensity_HbO = Conc_matrix(1:size(G,2),:);
intensity_HbR = Conc_matrix(size(G,2)+1:end,:);

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

first_stim = stim(1).data(1,:);
t = acquired.data.time;

figure
subplot(321)
time = first_stim - 1;
scatter3(brain_vertices_masked(:,1),brain_vertices_masked(:,2),brain_vertices_masked(:,3),1,intensity_HbO(:,find(t>time,1)))
colorbar
caxis([-1 1]*1e-3)
title('-1s')
view(60,7)

subplot(322)
time = first_stim + 0;
scatter3(brain_vertices_masked(:,1),brain_vertices_masked(:,2),brain_vertices_masked(:,3),1,intensity_HbO(:,find(t>time,1)))
colorbar
caxis([-1 1]*1e-3)
title('0s')
view(60,7)

subplot(323)
time = first_stim + 5;
scatter3(brain_vertices_masked(:,1),brain_vertices_masked(:,2),brain_vertices_masked(:,3),1,intensity_HbO(:,find(t>time,1)))
colorbar
caxis([-1 1]*1e-3)
title('5s')
view(60,7)

subplot(324)
time = first_stim + 10;
scatter3(brain_vertices_masked(:,1),brain_vertices_masked(:,2),brain_vertices_masked(:,3),1,intensity_HbO(:,find(t>time,1)))
colorbar
caxis([-1 1]*1e-3)
title('10s')
view(60,7)

subplot(325)
time = first_stim + 15;
scatter3(brain_vertices_masked(:,1),brain_vertices_masked(:,2),brain_vertices_masked(:,3),1,intensity_HbO(:,find(t>time,1)))
colorbar
caxis([-1 1]*1e-3)
title('15s')
view(60,7)

subplot(326)
time = first_stim + 20;
scatter3(brain_vertices_masked(:,1),brain_vertices_masked(:,2),brain_vertices_masked(:,3),1,intensity_HbO(:,find(t>time,1)))
colorbar
caxis([-1 1]*1e-3)
title('20s')
view(60,7)

% figure
% subplot(321)
% [h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>-1,1)), light_onoff, visible, axes_order);
% colorbar
% % caxis([-8 100])
% title('-1s')
% 
% subplot(322)
% [h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>0,1)), light_onoff, visible, axes_order);
% colorbar
% % caxis([-8 100])
% title('0s')
% 
% subplot(323)
% [h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>5,1)), light_onoff, visible, axes_order);
% colorbar
% % caxis([-8 100])
% title('5s')
% 
% subplot(324)
% [h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>10,1)), light_onoff, visible, axes_order);
% colorbar
% % caxis([-8 100])
% title('10s')
% 
% subplot(325)
% [h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>15,1)), light_onoff, visible, axes_order);
% colorbar
% % caxis([-8 100])
% title('15s')
% 
% subplot(326)
% [h, l] = displayIntensityOnMesh(mesh, intensity_HbO(:,find(tHRF>20,1)), light_onoff, visible, axes_order);
% colorbar
% % caxis([-8 100])
% title('20s')
% 
% set(gcf,'position',[ 1001         135        1147        1204])