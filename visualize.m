function visualize(Conc,acquired_path,M,Proc_data)

At_file = 'atlasViewer.mat';
atlasViewer = load([acquired_path,At_file]);
warning off

axes_order = [2,1,3];
faces = atlasViewer.fwmodel.mesh.faces;
brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
intensity = zeros(size(brain_vertices,1), size(Conc.intensity_HbO,2));
intensity(M.mask_brain,:) = Conc.intensity_HbO;
% 
% first_stim = Proc_data.stim(1).data(1,1);
t = Proc_data.tHRF;

figure
% subplot(321)
% t_in_tHRF = -1;
% time = t_in_tHRF;
% position = find(t>time,1);
% trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
%           intensity(:,position),'facecolor','interp','edgealpha',0, 'visible','on');
% colorbar
% caxis([-1 1]*1e-7)
% title('-1s')
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([128.0, 2238.8, 130.0])
% camup([-1.0, 0.0, 0.0])

% subplot(322)
% t_in_tHRF = 0;
% time = t_in_tHRF;
% position = find(t>time,1);
% trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
%           intensity(:,position),'facecolor','interp','edgealpha',0, 'visible','on');
% colorbar
% caxis([-1 1]*1e-7)
% title('0s')
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([128.0, 2238.8, 130.0])
% camup([-1.0, 0.0, 0.0])

% subplot(323)
% t_in_tHRF = 5;
% time = t_in_tHRF;
% position = find(t>time,1);
% trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
%           intensity(:,position),'facecolor','interp','edgealpha',0, 'visible','on');
% colorbar
% caxis([-1 1]*1e-7)
% title('5s')
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([128.0, 2238.8, 130.0])
% camup([-1.0, 0.0, 0.0])

% subplot(324)
t_in_tHRF = 10;
time = t_in_tHRF;
position = find(t>time,1);
h = trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
          intensity(:,position),'facecolor','interp','edgealpha',0, 'visible','on');
set(h,'diffusestrength',.9,'specularstrength',.12,'ambientstrength',.2);
colorbar
caxis([-0.5 1]*0.8e-5)
title('10s')
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])

if(~exist('light_onoff') | (exist('light_onoff') & strcmp(light_onoff,'on')))
    l = camlight;
    set(l,'Position',[50 2000 100]);

    l2 = camlight;
    set(l2,'Position',[50 -100 -100]);

    camlight(0,0);
end
lighting phong;
colormap jet
% subplot(325)
% t_in_tHRF = 15;
% time = t_in_tHRF;
% position = find(t>time,1);
% trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
%           intensity(:,position),'facecolor','interp','edgealpha',0, 'visible','on');
% colorbar
% caxis([-1 1]*1e-7)
% title('15s')
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([128.0, 2238.8, 130.0])
% camup([-1.0, 0.0, 0.0])
% 
% subplot(326)
% t_in_tHRF = 20;
% time = t_in_tHRF;
% position = find(t>time,1);
% trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
%           intensity(:,position),'facecolor','interp','edgealpha',0, 'visible','on');
% colorbar
% caxis([-1 1]*1e-7)
% title('20s')
% view(-90,0)
% camtarget([128.0, 132.0, 130.0])
% campos([128.0, 2238.8, 130.0])
% camup([-1.0, 0.0, 0.0])

% set(gcf,'position',[ 10         10        1468        1204])