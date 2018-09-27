function L = logMap(p,x)
%     [u,G] = eig(p);
%     g = u*sqrt(G);
%     gi = diag(1./sqrt(diag(G))) * u';
%     y = gi * x * gi';
%     [v,S] = eig(y);
%     gv = g * v;
%     L = gv * diag(log(diag(S))) * gv';

ph = p^0.5;
phi = p^-0.5;

tmp = phi*x*phi;
tmp = (tmp + tmp') / 2;
[u,G] = eig(tmp);

tmp = u*diag(log(diag(G)))*u';
tmp = (tmp + tmp') / 2;

L = ph*tmp*ph;
L = (L + L') / 2;

end