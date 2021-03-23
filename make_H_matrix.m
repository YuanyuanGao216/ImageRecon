function [HTH,HTY] = make_H_matrix(A, G, T, Y, device, n_batch_brain, n_batch_scalp)

canUseGPU=true;

if strcmp(device, 'gpu')
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
%     if ~isfile('H_brain.mat')
    if 1
        fprintf('Calculating brain\n')
        n_batch = n_batch_brain;
        H_finished = false;
        i = 1;
        while i <= 100
            try
                fprintf('number of batch is %d\n',n_batch)
                H_brain = calculate_H_batch_gpu(A.Amatrix_brain, G.G_brain, T.T_HbO_brain, T.T_HbR_brain, Y, n_batch);
                fprintf('finished H %d\n',n_batch)
                H_finished = true;
                save('H_brain.mat','H_brain')
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
    else
        load('H_brain.mat','H_brain')
    end
%     if ~isfile('H_scalp.mat')
    if 1
        fprintf('Calculating scalp\n')
        n_batch = n_batch_scalp;
        H_finished = false;
        i = 1;
        while i <= 100
            try
                fprintf('number of batch is %d\n',n_batch)
                H_scalp = calculate_H_batch_gpu(A.Amatrix_scalp, G.G_scalp, T.T_HbO_scalp, T.T_HbR_scalp, Y, n_batch);
                fprintf('finished H %d\n',n_batch)
                H_finished = true;
                save('H_scalp.mat','H_scalp')
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
    else
        load('H_scalp.mat','H_scalp')
    end
    
    Y = gpuArray(Y);
    HTH = zeros(size(H_brain,2));
    HTY = zeros(size(H_scalp,2),1);
    chunk = 100;
    for i = 1:chunk:size(H_brain,2)
        index_start = i;
        index_end = i + chunk - 1;
        if index_end > size(H_brain,2)
            index_end = size(H_brain,2);
        end
        index1 = index_start : index_end;
        H_chunk_1 = gpuArray(H_brain(:,index1));
        
        hty = H_chunk_1'*Y;
        HTY(index1,:) = hty;
        for j = 1:chunk:size(H_brain,2)
            index_start = j;
            index_end = j + chunk - 1;
            if index_end > size(H_brain,2)
                index_end = size(H_brain,2);
            end
            index2 = index_start : index_end;
            H_chunk_2 = gpuArray(H_brain(:,index2));
            
            hth = H_chunk_1'*H_chunk_2;
            HTH(index1,index2) = hth;
        end
        for j = 1:chunk:size(H_scalp,2)
            index_start = j;
            index_end = j + chunk - 1;
            if index_end > size(H_scalp,2)
                index_end = size(H_scalp,2);
            end
            index2 = index_start : index_end;
            H_chunk_2 = gpuArray(H_scalp(:,index2));
            
            hth = H_chunk_1'*H_chunk_2;
            HTH(index1,size(H_brain,2) + index2) = hth;
            HTH(size(H_brain,2) + index2,index1) = hth';
        end
        
    end
    for i = 1:chunk:size(H_scalp,2)
        index_start = i;
        index_end = i + chunk - 1;
        if index_end > size(H_scalp,2)
            index_end = size(H_scalp,2);
        end
        index1 = index_start : index_end;
        H_chunk_1 = gpuArray(H_scalp(:,index1));
        
        hty = H_chunk_1'*Y;
        HTY(size(H_brain,2) + index1,:) = hty;
        for j = 1:chunk:size(H_scalp,2)
            index_start = j;
            index_end = j + chunk - 1;
            if index_end > size(H_scalp,2)
                index_end = size(H_scalp,2);
            end
            index2 = index_start : index_end;
            H_chunk_2 = gpuArray(H_scalp(:,index2));
            
            hth = H_chunk_1'*H_chunk_2;
            HTH(size(H_brain,2) + index1 , size(H_brain,2) + index2) = hth;
        end
    end
%     H = gpuArray(H);
%     Y = gpuArray(Y);
%     HTY = (H'*Y);
%     HTH = H'*H;
end

if strcmp(device, 'cpu') || ~canUseGPU
    [HTH,HTY] = calculate_H_batch_cpu(A, G, T, Y);
end
end

function H = calculate_H_batch_gpu(A, G, T_HbO, T_HbR, Y, n_batch)

n_g = size(G,1);
n_t = size(T_HbO,2);

n_g_batch = floor(n_g/n_batch);
if n_g_batch == 0
    n_g_batch = 1;
end
H = zeros(length(Y),n_g*n_t);

Amatrix = gpuArray(A);

% f = waitbar(0,'1','Name','Calculating H for brain',...
%     'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
% 
% setappdata(f,'canceling',0);
i_g = 1;
while i_g <= n_g
%     if getappdata(f,'canceling')
%         break
%     end
%     waitbar(i_g/n_g,f,sprintf('%0.3f%%',i_g/n_g*100))
    
    g_start     =   i_g;
    g_end       =   i_g + n_g_batch - 1 ;
    if g_end > n_g
        g_end = n_g;
    end
    g_index     =   g_start:g_end;
    
    G_brain = gpuArray(G(g_index,:));
    T_HbO = gpuArray(T_HbO);
    T_HbR = gpuArray(T_HbR);
    
    HbO = kron(reshape(G_brain',[],1),reshape(T_HbO,1,[]));
    HbR = kron(reshape(G_brain',[],1),reshape(T_HbR,1,[]));
    [m, n] = size(G_brain);
    [k, l] = size(T_HbO);
    for j = 1:m
        for i = 1:l
            HbO_sub = HbO((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
            HbR_sub = HbR((j-1)*n + 1 : j*n, (i-1)*k + 1: i*k);
            gt = [HbO_sub;HbR_sub];
            Agt = Amatrix*gt;
            Agt_cpu = gather(Agt);
            H(:,(g_start + j- 2) * l + i ) = reshape(Agt_cpu,[],1);
        end
        
    end
    
    i_g = i_g + n_g_batch;
end
% delete(f)
end

function [HTH,HTY] = calculate_H_batch_cpu(A, G, T, Y)

n_g = size(G.G_brain,1);
n_t = size(T.T_HbO_brain,2);
H_brain = zeros(length(Y),n_g*n_t);

f = waitbar(0,'1','Name','Calculating H for brain',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

j = 1;
for i_g = 1:n_g
    if getappdata(f,'canceling')
        break
    end
    waitbar(i_g/n_g,f,sprintf('%0.3f%%',i_g/n_g*100))
    for i_t = 1:n_t
        HbO = kron(G.G_brain(i_g,:),T.T_HbO_brain(:,i_t));
        HbR = kron(G.G_brain(i_g,:),T.T_HbR_brain(:,i_t));
        gt = [HbO,HbR]';
        Agt = A.Amatrix_brain*gt;
        H_brain(:,j) = reshape(Agt,[],1);
        j = j + 1;
    end
end
delete(f)

n_g = size(G.G_scalp,1);
n_t = size(T.T_scalp,2);
H_scalp = zeros(length(Y),n_g*n_t);

f = waitbar(0,'1','Name','Calculating H for brain',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');

setappdata(f,'canceling',0);

j = 1;
for i_g = 1:n_g
    if getappdata(f,'canceling')
        break
    end
    waitbar(i_g/n_g,f,sprintf('%0.3f%%',i_g/n_g*100))
    for i_t = 1:n_t
        HbO = kron(G.G_scalp(i_g,:),T.T_scalp(:,i_t));
        gt = [HbO,HbO]';
        Agt = A.Amatrix_scalp*gt;
        H_scalp(:,j) = reshape(Agt,[],1);
        j = j + 1;
    end
end
delete(f)
H = [H_brain,H_scalp];
save('H_CPU.mat','H')
HTY = (H'*Y);
HTH = H'*H;

end