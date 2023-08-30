function [lonOut,latOut]=startangle2end(lon,lat,yaw,s)
%根据p1以及航向和距离，计算终点p2的位置以及方位角
%输入：lon,lat 点1经纬度 ，yaw点的方位角 deg  s 点沿方位角的距离
%输出：终点的经纬度
lon = vpa(lon,10);
lat = vpa(lat,10);
lon1 = lon*pi/180 ; 
lat1 = lat*pi/180;  %p1	
alpha_1 = deg2rad(yaw + 360); %  方向角

f = 1 /  298.257223563;
a= 6378137.0;	
b= 6356752.314245;

sin_alpha1 = sin(alpha_1);
cos_alpha1 = cos(alpha_1);

tanU1 = (1-f) * tan(lat1); cosU1 = 1 / sqrt((1 + tanU1*tanU1)); sinU1 = tanU1 * cosU1;
delta_1 = atan2(tanU1, cos_alpha1);
sin_alpha = cosU1 * sin_alpha1;
cosSq_alpha = 1 - sin_alpha*sin_alpha;
uSq = cosSq_alpha * (a*a - b*b) / (b*b);
A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));

delta = s / (b*A);
dleta_ = 0;
while abs(delta-dleta_) > 1e-12
    cos2deltaM = cos(2*delta_1 + delta);
    sin_delta = sin(delta);
    cos_delta = cos(delta);
    delta_delta = B*sin_delta*(cos2deltaM+B/4*(cos_delta*(-1+2*cos2deltaM*cos2deltaM)-B/6*cos2deltaM*(-3+4*sin_delta*sin_delta)*(-3+4*cos2deltaM*cos2deltaM)));
    dleta_ = delta;
    delta = s / (b*A) + delta_delta;
end

tmp = sinU1*sin_delta - cosU1*cos_delta*cos_alpha1;
lat2 = atan2(sinU1*cos_delta + cosU1*sin_delta*cos_alpha1, (1-f)*sqrt(sin_alpha*sin_alpha + tmp*tmp)); %目标点纬度
lon =  atan2(sin_delta*sin_alpha1, cosU1*cos_delta - sinU1*sin_delta*cos_alpha1);
C = f/16*cosSq_alpha*(4+f*(4-3*cosSq_alpha));
L = lon - (1-C) * f * sin_alpha *(delta + C*sin_delta*(cos2deltaM+C*cos_delta*(-1+2*cos2deltaM*cos2deltaM)));
lon2 =rem((lon1+L+3*pi) ,(2*pi)) - pi; 
% lon2 =roundn(rem((lon1+L+3*pi) ,(2*pi)) - pi, -10); %normalise to -180...+180 目标点精度

lonOut = vpa(lon2*180/pi);
latOut = vpa(lat2*180/pi);

revAz = atan2(sin_alpha, -tmp); %最终方位角

end