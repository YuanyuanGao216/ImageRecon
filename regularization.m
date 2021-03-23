function [Adot, Adot_scalp] = regularization(Adot, Adot_scalp)

% spatial regularization on brain only
for j = 1:size(Adot,3)
    J = single(squeeze(Adot(:,:,j)));
    try
        JTJ = diag(J'*J);
    catch
        close(h);
        menu(sprintf('Out of memory: JTJ = diag(J''*J) generates matrix that is too large!!'), 'Okay');
        return;
    end
    L = beta.*max(JTJ);
    LL = sqrt(JTJ+L)';
    LL = diag(1./LL);
    % Normalise J
    Adot(:,:,j) = J*LL;
end

[u1,s1,v1]=svds(double([squeeze(Adot(:,:,1)) squeeze(Adot_scalp(:,:,1))]),size(Adot,1)); max_sing1 = max(s1(:));
[u2,s2,v2]=svds(double([squeeze(Adot(:,:,2)) squeeze(Adot_scalp(:,:,2))]),size(Adot,1)); max_sing2 = max(s2(:));

% regularization parameters
alpha1 = alpha * max_sing1 ;
alpha2 = alpha * max_sing2 ;

% new sensitivity with regularization
Anew(:,:,1) = u1 * sqrtm(s1*s1 + alpha1^2*eye(size(s1))) *v1';
Anew(:,:,2) = u2 * sqrtm(s2*s2 + alpha2^2*eye(size(s2))) *v2';