function x = snapzero(x, tol)
    if nargin<2 
        tol = 1e-12; 
    end
    x(abs(x) < tol) = 0;
end
