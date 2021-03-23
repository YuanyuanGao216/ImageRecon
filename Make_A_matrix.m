function A = Make_A_matrix(Adot, Adot_scalp, E, M, alpha, beta)

mask = M.mask_brain;
mask_scalp = M.mask_scalp;

Adot = Adot(:, mask, :);
Adot_scalp = Adot_scalp(:, mask_scalp, :);

% [Adot, Adot_scalp] = regularization(Adot, Adot_scalp);

    
Amatrix_brain = [squeeze(Adot(:,:,1))*E(1,1) squeeze(Adot(:,:,1))*E(1,2);
            squeeze(Adot(:,:,2))*E(2,1) squeeze(Adot(:,:,2))*E(2,2)]; % 100, 40008
Amatrix_scalp = [squeeze(Adot_scalp(:,:,1))*E(1,1) squeeze(Adot_scalp(:,:,1))*E(1,2);
            squeeze(Adot_scalp(:,:,2))*E(2,1) squeeze(Adot_scalp(:,:,2))*E(2,2)]; % 100, 19124
        
A.Amatrix_brain = Amatrix_brain;
A.Amatrix_scalp = Amatrix_scalp;