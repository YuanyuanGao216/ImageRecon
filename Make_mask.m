function M = Make_mask(At_file,threshold, Adot, Adot_scalp)

atlasViewer = load(At_file);

brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
scalp_vertices = atlasViewer.fwmodel.mesh_scalp.vertices; %9563*3

%% making mask here, output is mask and mask_scalp
figure('name','brain sensitivity')
axes_order = [2,1,3];
intensity = log10(sum(Adot(:,:,1),1));
faces = atlasViewer.fwmodel.mesh.faces;
trisurf(faces, brain_vertices(:,axes_order(1)), brain_vertices(:,axes_order(2)), brain_vertices(:,axes_order(3)), ...
          intensity,'facecolor','interp','edgealpha',0, 'visible','on');
colorbar
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
mask = find(intensity>threshold);
fprintf('length is %d\n', length(mask))
set(gcf,'position',[10 10 560 420])

figure('name','scalp sensitivity')
intensity_scalp = log10(sum(Adot_scalp(:,:,1),1));
faces_scalp = atlasViewer.fwmodel.mesh_scalp.faces;
trisurf(faces_scalp, scalp_vertices(:,axes_order(1)), scalp_vertices(:,axes_order(2)), scalp_vertices(:,axes_order(3)), ...
          intensity_scalp,'facecolor','interp','edgealpha',0, 'visible','on');
% trisurf(faces_scalp,scalp_vertices(:,1),scalp_vertices(:,2),scalp_vertices(:,3),intensity_scalp,'edgecolor','none')
colorbar
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
mask_scalp = find(intensity_scalp>threshold);
fprintf('scalp length is %d\n', length(mask_scalp))
set(gcf,'position',[10 10 560 420])

% here i want to show what is left after mask
% figure
% trisurf(T,brain_vertices(:,1),brain_vertices(:,2),brain_vertices(:,3),'edgecolor',[0.1 0.1 0.1],'edgecolor','none')
% hold on
% trisurf(T(),brain_vertices(mask,1),brain_vertices(mask,2),brain_vertices(mask,3),'edgecolor',[0.1 0.1 0.1],'edgecolor','none')

M.mask_brain = mask;
M.mask_scalp = mask_scalp;


