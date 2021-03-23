function G = Make_G_matrix(At_file,M, threshold_brain, threshold_scalp, sigma_brain, sigma_scalp)

atlasViewer = load(At_file);

brain_vertices = atlasViewer.fwmodel.mesh.vertices; %20004*3
scalp_vertices = atlasViewer.fwmodel.mesh_scalp.vertices; %9563*3

mask_brain = M.mask_brain;
mask_scalp = M.mask_scalp;

% determine kernel positions
brain_vertices_masked = brain_vertices(mask_brain,:);
scalp_vertices_masked = scalp_vertices(mask_scalp,:);


[brain_vertices_new] = down_sample_vertices(brain_vertices_masked, threshold_brain);
[scalp_vertices_new] = down_sample_vertices(scalp_vertices_masked, threshold_scalp);

% visualize the new vertices on old mesh, trimesh
axes_order = [2,1,3];

figure('name','brain')
faces = atlasViewer.fwmodel.mesh.faces;
trimesh(faces,brain_vertices(:,axes_order(1)),brain_vertices(:,axes_order(2)),brain_vertices(:,axes_order(3)))
hold on
plot3(brain_vertices_new(:,axes_order(1)),brain_vertices_new(:,axes_order(2)),brain_vertices_new(:,axes_order(3)),'ro');
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
set(gcf,'position',[10 10 560 420])

figure('name','scalp')
faces_scalp = atlasViewer.fwmodel.mesh_scalp.faces;
trimesh(faces_scalp, scalp_vertices(:,axes_order(1)),scalp_vertices(:,axes_order(2)),scalp_vertices(:,axes_order(3)))
hold on
plot3(scalp_vertices_new(:,axes_order(1)),scalp_vertices_new(:,axes_order(2)),scalp_vertices_new(:,axes_order(3)),'ro');
view(-90,0)
camtarget([128.0, 132.0, 130.0])
campos([128.0, 2238.8, 130.0])
camup([-1.0, 0.0, 0.0])
set(gcf,'position',[10 10 560 420])

% make kernels

G_brain = make_kernel_matrix_gpu(brain_vertices_new, brain_vertices_masked,sigma_brain);
G_scalp = make_kernel_matrix_gpu(scalp_vertices_new, scalp_vertices_masked,sigma_scalp);

G.G_brain = G_brain;
G.G_scalp = G_scalp;