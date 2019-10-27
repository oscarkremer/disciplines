%%%% faz o plot das clusterizacoess %%%%

function plotClass(X, label, iter)
[d,n] = size(X);
if nargin == 1
    label = ones(n,1);
end
assert(n == length(label));

color = 'brgmcyk';
m = length(color);
c = max(label);

figure(iter);
clf;
hold on;
switch d
    case 2
        view(2);
        for i = 1:c
            idc = label==i;
            scatter(X(1,idc),X(2,idc),10,color(mod(i-1,m)+1));
        end
    case 3
        view(3);
        for i = 1:c
            idc = label==i;
            scatter3(X(1,idc),X(2,idc),X(3,idc),10,color(mod(i-1,m)+1));
        end
    otherwise
        error('ERROR');
end
axis equal;
grid on;
hold off;