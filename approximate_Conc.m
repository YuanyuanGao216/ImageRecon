function Conc = approximate_Conc(G,T,beta, device, n_batch_brain)
canUseGPU=true;

if strcmp(device,'gpu')
    try
        fprintf('Finding GPU...\n')
        if gpuDeviceCount > 1
            for ii = 1:gpuDeviceCount
                g = gpuDevice(ii);
                fprintf(1,'Device %i has ComputeCapability %s \n', ...
                        g.Index,g.ComputeCapability)
            end
            prompt = 'Which GPU? ';
            x = input(prompt);
            d = gpuDevice(x);
        else
            d = gpuDevice;
        end
        display(d)
    catch
        canUseGPU=false;
        fprintf('no GPU. Using cpu now\n')
    end
    fprintf('Calculating brain\n')
    n_batch = n_batch_brain;
    H_finished = false;
    i = 1;
    while i <= 100
        try
            fprintf('number of batch is %d\n',n_batch)
            Conc = calculate_Conc(G.G_brain, T.t_HbO_brain, T.t_HbR_brain, beta, n_batch);
            fprintf('finished Conc %d\n',n_batch)
            H_finished = true;
%             save('Conc.mat','Conc')
            break
        catch
            n_batch = n_batch*2;
            i = i + 1;
            continue
        end

    end
    if H_finished == false
        fprintf('Cannot fit into GPU memory even with %d batches\n', n_batch)
        return
    end
    
elseif strcmp(device, 'cpu') || ~canUseGPU
    n_g = size(G.G_brain,1);
    n_t = size(T.t_HbO_brain,2);
    beta_brain = beta(1:n_g*n_t);

    Conc_matrix = zeros(size(G.G_brain,2)*2,size(T.t_HbO_brain,1));
    G_brain = G.G_brain;
    T_HbO_brain = T.t_HbO_brain;
    T_HbR_brain = T.t_HbR_brain;

    j = 1;

    for i_g = 1:n_g
        for i_t = 1:n_t
            HbO = kron(G_brain(i_g,:),T_HbO_brain(:,i_t));
            HbR = kron(G_brain(i_g,:),T_HbR_brain(:,i_t));
            gt = [HbO,HbR]';
            Conc_matrix = Conc_matrix + gt.*beta_brain(j);
            j = j + 1;
        end
    end

    intensity_HbO = Conc_matrix(1:size(G.G_brain,2),:);
    intensity_HbR = Conc_matrix(size(G.G_brain,2)+1:end,:);

    Conc.intensity_HbO = intensity_HbO;
    Conc.intensity_HbR = intensity_HbR;
end
end

function Conc = calculate_Conc(G, T_HbO_brain, T_HbR_brain, beta, n_batch)
n_g = size(G,1);
n_t = size(T_HbO_brain,2);
beta_brain = beta(1:n_g*n_t);

Conc_matrix = zeros(size(G,2)*2,size(T_HbO_brain,1),'gpuArray');

n_g_batch = floor(n_g/n_batch);

i_g = 1;
while i_g <= n_g
    g_start     =   i_g;
    g_end       =   i_g + n_g_batch - 1 ;
    if g_end > n_g
        g_end = n_g;
    end
    g_index     =   g_start:g_end;
    
    G_brain = gpuArray(G(g_index,:));
    T_HbO = gpuArray(T_HbO_brain);
    T_HbR = gpuArray(T_HbR_brain);
    beta = beta_brain(g_index);
    
    HbO = kron(reshape(G_brain',[],1),reshape(T_HbO,1,[]));
    HbR = kron(reshape(G_brain',[],1),reshape(T_HbR,1,[]));
    [m, n] = size(G_brain);
    [k, l] = size(T_HbO);
    h = 1;
    for j = 1:m
        for i = 1:l
            HbO_sub = HbO((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
            HbR_sub = HbR((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
            gt = [HbO_sub;HbR_sub];
            Conc_matrix = Conc_matrix + gt.*beta(h);
            h = h + 1;
        end
    end
    i_g = i_g + n_g_batch;
end

intensity_HbO = Conc_matrix(1:size(G,2),:);
intensity_HbR = Conc_matrix(size(G,2)+1:end,:);

Conc.intensity_HbO = intensity_HbO;
Conc.intensity_HbR = intensity_HbR;

Conc = gather(Conc);

end