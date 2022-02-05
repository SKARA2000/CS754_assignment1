function [out_f] = omp_reconstruction(S, y, eps, patch_size)
    % Function reconstructs the image using the coded aperture and the measured 
    % coded snapshot. 
    % S : coded Aperture Values
    % y : Coded snapshot
    % eps : error bound for the OMP reconstructin algorithm
    % patch_size : Size of the patches to be considered for reconstruction
    % out_f : reconstructed image
    T = size(S, 3);
    DCT_basis = dctmtx(patch_size);
    DCT_2d_mtx = kron(DCT_basis, DCT_basis);
    psi = kron(eye(T), DCT_2d_mtx);

    out_f = zeros(size(S, 1), size(S, 2), size(S, 3));
    overlaps = zeros(size(S, 1), size(S, 2), size(S, 3));
    
    count = 1;
    num_patches = (size(y, 1) - patch_size + 1)*(size(y, 2) - patch_size + 1);
    for i = 1:(size(y, 1) - patch_size + 1)
        for j = 1:(size(y, 2) - patch_size + 1)
            S_ij = S(i:(i + patch_size - 1), j:(j + patch_size - 1), :);
            y_ij = y(i:(i + patch_size - 1), j:(j + patch_size - 1));
            phi = zeros(patch_size^2, (patch_size^2)*T);
            % Rearrange the coded aperture values to a block-diagonal matrix
            % with its diagonals being the coded aperture used for each frame
            for t = 1:T
                S_ijt = S_ij(:, :, t);
                phi(:, ((t-1)*(patch_size^2) + 1):(t*(patch_size^2))) = diag(reshape(S_ijt, patch_size^2, 1));
            end
            A_ij = phi * psi';
            b_ij = reshape(y_ij, (patch_size^2), 1);
            % Perform the OMP algorithm with the error bound
            [theta_ij] = omp(A_ij, b_ij, eps);

            reconstructed_x_ij = psi' * theta_ij;
            for t = 1:T
                temp_t = reshape(reconstructed_x_ij(((t - 1)*(patch_size^2) + 1): (t*(patch_size^2)), :), patch_size, patch_size);
                out_f(i:(i + patch_size - 1), j:(j + patch_size - 1), t) = out_f(i:(i + patch_size - 1), j:(j + patch_size - 1), t) + temp_t;
                overlaps(i:(i + patch_size - 1), j:(j + patch_size - 1), t) = overlaps(i:(i + patch_size - 1), j:(j + patch_size - 1), t) + ones(patch_size, patch_size);
            end
            count = count + 1;
            if rem(count, 1000) == 0
                fprintf("Done with (%d/%d) patches\n", count, num_patches);
            end
        end
    end
    fprintf("Done with all patches\n");
    out_f = out_f./overlaps;
end 

function [x] = omp(A, b, eps)
    r = b;
    x = zeros(size(A, 2), 1);
    T = [];
    i = 0;
    while((norm(r) > eps) && (i < 64))
        % Selecting the column with the maximal dot product with r
        nA = sqrt(sum(A.*A, 1));
        numerator = sum(A.*repmat(r, 1, size(A, 2)), 1);
        [~, j] = max(abs(numerator./nA));
        % Adding to the support set
        T = [T j];
        A_T = A(:, T);
        i = i + 1;
        % Update x and the residual vector
        x(T) = pinv(A_T)*b;
        r = b - A_T*x(T);
    end
end