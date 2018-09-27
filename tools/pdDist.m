function d = pdDist(x,y)
   t = x^-0.5*y*x^-0.5;
   [u,G] = eig(t);
   d = sum(log(diag(G)).^2);
end