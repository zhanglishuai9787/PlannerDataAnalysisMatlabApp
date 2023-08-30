function d = SphereDist(x,y,R)
%x为A点[经度, 纬度], y为B点[经度, 纬度]
if nargin < 3
    R = 6378.137;
end
x = x*pi/180;
y = y*pi/180;
DeltaS = acos(cos(x(2))*cos(y(2))*cos(x(1)-y(1))+sin(x(2))*sin(y(2)));
d = R*DeltaS;

end
