function [dx, dy] = subpixel_quadratic_fit(patch)
    % 3x3 quadratic fitting
    [X, Y] = meshgrid(-1:1, -1:1);
    Z = patch(:);
    A = [X(:).^2, Y(:).^2, X(:).*Y(:), X(:), Y(:), ones(9,1)];
    coeff = A\Z;
    a=coeff(1); b=coeff(2); c=coeff(3); d=coeff(4); e=coeff(5);
    H = [2*a, c; c, 2*b];
    g = [d; e];
    offset = -H\g;
    dx = offset(1);
    dy = offset(2);
end