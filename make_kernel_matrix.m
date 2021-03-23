function Kernel_matrix = make_kernel_matrix(brain_vertices_new, brain_vertices,sigma)

Cov_matrix = sigma^2 * eye(3);
n_kernels = size(brain_vertices_new, 1);
n_vertices = size(brain_vertices, 1);
Kernel_matrix = zeros(n_kernels, n_vertices);

for i=1:n_kernels
    mu = brain_vertices_new(i,:);
    kernel_vector = zeros(1,n_vertices);
    for j = 1:n_vertices
        x = brain_vertices(j,:);
        kernel_vector(j) = exp(-(x-mu)/Cov_matrix*(x-mu)'/2)/(sqrt((2*pi)^3*det(Cov_matrix)));
    end
    kernel_vector = kernel_vector./sum(kernel_vector);
    Kernel_matrix(i,:) = kernel_vector;
end
