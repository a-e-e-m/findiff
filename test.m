    dim = 35;

    a = ones(dim,1);
    b = ones(dim,1);
    c = 4 * ones(dim,1);
    d = ones(dim,1);
    e = ones(dim,1);

    Diag = [a,b,c,d,e];

    d = [-6;-1;0;1;6];

    A = spdiags(Diag,d,15,5);
    asa=full(A)