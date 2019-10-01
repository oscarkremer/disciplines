function J = distortion(mu1, mu2, class1, class2)
  m1 = size(class1); m1 = m1(1);
  m2 = size(class2); m2 = m2(1);
  dist1 = 0;
  dist2 = 0;
  for i=1:m1
    dist1 = dist1 + (class1(i,1) - mu1(1))^2 + (class1(i,2) - mu1(2))^2;
  end
  for i=1:m2
    dist2 = dist2 + (class2(i,1) - mu2(1))^2 + (class2(i,2) - mu2(2))^2;
  end
  J = dist1 + dist2;
end