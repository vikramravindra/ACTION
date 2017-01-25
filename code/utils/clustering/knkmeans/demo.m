clear; close all;
d = 2;
k = 3;
n = 500;
[X,label] = kmeansRnd(d,k,n);
init = ceil(k*rand(1,n));
[y,mse,model] = knKmeans(X,init,@knLin);
plotClass(X,y)

idx = 1:2:n;
Xt = X(:,idx);
t = knKmeansPred(model, Xt);
plotClass(Xt,t)