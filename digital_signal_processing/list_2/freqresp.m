function [H ] = freqresp(b, a, w)
  num_b = size(b); num_b = num_b(1);
  dem_a = size(a); dem_a = dem_a(1);
  discrete_w = size(w); discrete_w = discrete_w(1);
  H = zeros(discrete_w, 1);
  for j=1:discrete_w
    num = 1;
    dem = 1;
    for m=1:num_b
      num = num*b(m)*exp(-i*w(j)*(m-1));
    end
    for l=1:dem_a
      dem = dem*a(l)*exp(-i*w(j)*(l-1));
    H(j) = num/dem;
    end
  end
 end