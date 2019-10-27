%%%%  %%%%

function s = logsumexp(X, dim)

  if nargin == 1, 
    dim = find(size(X)~=1,1);
    if isempty(dim), dim = 1; end
  end

  y = max(X,[],dim);
  s = y+log(sum(exp(bsxfun(@minus,X,y)),dim));  
  i = isinf(y);
  if any(i(:))
      s(i) = y(i);
  end