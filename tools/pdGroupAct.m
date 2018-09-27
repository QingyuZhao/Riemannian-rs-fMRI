function Y = pdGroupAct(X,p1,p2)
    [u,G] = eig(p1);
    gi = diag(1./sqrt(diag(G))) * u';
    
    [u,G] = eig(p2);
    f = u*sqrt(G);
    
    fgi = f*gi;
    Y = fgi*X*fgi';
end