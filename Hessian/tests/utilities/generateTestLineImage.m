function I = generateTestLineImage(sz,angleDeg,width)
[x,y] = meshgrid(1:sz,1:sz);
cx = sz/2; cy = sz/2;
theta = deg2rad(angleDeg);
xr = (x-cx)*cos(theta)+(y-cy)*sin(theta);
I = zeros(sz);
I(abs(xr)<=width)=1;
I = imgaussfilt(I,1);
end
