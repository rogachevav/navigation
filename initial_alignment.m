function [L] = initial_alignment(A, phi, h, a, e2, u)
L = zeros(3, 3);
means = mean(A(:,2:end));
L(:, 3) = means(4:end) / g_phi(phi, h, a, e2);
L(:, 2) = (means(1:3) - (L(:, 3)' * u * sin(phi))) / (u * cos(phi));
L(:, 1) = cross(L(:, 2), L(:, 3))';
L(:, 3) =L(:, 3)/ norm(L(:, 3), 2);
L(:, 2) =L(:, 2)/ norm(L(:, 2), 2);
L(:, 1) =L(:, 1)/ norm(L(:, 1), 2);
end