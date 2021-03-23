function Kernel_matrix = make_kernel_matrix_gpu(brain_vertices_new, brain_vertices,sigma)


try
    fprintf('finding GPU for G matrix\n')
    d = gpuDevice;
    display(d)
catch
    fprintf('no GPU. Using cpu now for G\n')
    Kernel_matrix = make_kernel_matrix(brain_vertices_new, brain_vertices_masked,sigma_brain);
end

brain_vertices_new = gpuArray(brain_vertices_new);
brain_vertices = gpuArray(brain_vertices);

Cov_matrix = gpuArray(sigma^2 * eye(3));
n_kernels = size(brain_vertices_new, 1);
n_vertices = size(brain_vertices, 1);
Kernel_matrix = zeros(n_kernels, n_vertices);

for i=1:n_kernels
    mu = brain_vertices_new(i,:);
    
    x_minus_mu = brain_vertices - repmat(mu,n_vertices,1);
    kernel_vector = exp(-diag(x_minus_mu/Cov_matrix*x_minus_mu')/2)/(sqrt((2*pi)^3*det(Cov_matrix)));
%     for j = 1:n_vertices
%         x = brain_vertices(j,:);
%         kernel_vector(j) = exp(-(x-mu)/Cov_matrix*(x-mu)'/2)/(sqrt((2*pi)^3*det(Cov_matrix)));
%     end
    kernel_vector = kernel_vector./sum(kernel_vector);
    Kernel_matrix(i,:) = gather(kernel_vector);
end
