function [lonOut,latOut]=startangle2endhigh(lon,lat,yaw,s)
%根据p1以及航向和距离，计算终点p2的位置以及方位角
%输入：lon,lat 点1经纬度 ，yaw点的方位角 deg  s 点沿方位角的距离
%输出：终点的经纬度
lon = vpa(lon,10);
lat = vpa(lat,10);
lon1 = vpa(lon*pi/180); 
lat1 = vpa(lat*pi/180);  %p1	
alpha_1 = vpa(deg2rad(yaw + 360)); %  方向角

f = 1 /  298.257223563;
a= 6378137.0;	
b= 6356752.314245;

sin_alpha1 = vpa(sin(alpha_1));
cos_alpha1 = vpa(cos(alpha_1));

tanU1 = vpa((1-f) * tan(lat1)); cosU1 = vpa(1 / sqrt((1 + tanU1*tanU1))); sinU1 = vpa(tanU1 * cosU1);
delta_1 = vpa(atan2(tanU1, cos_alpha1));
sin_alpha = vpa(cosU1 * sin_alpha1);
cosSq_alpha = vpa(1 - sin_alpha*sin_alpha);
uSq = vpa(cosSq_alpha * (a*a - b*b) / (b*b));
A = vpa(1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq))));
B = vpa(uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq))));

delta = vpa(s / (b*A));
dleta_ = 0;
while vpa(abs(delta-dleta_)) > 1e-12
    cos2deltaM = vpa(cos(2*delta_1 + delta));
    sin_delta = vpa(sin(delta));
    cos_delta = vpa(cos(delta));
    delta_delta = vpa(B*sin_delta*(cos2deltaM+B/4*(cos_delta*(-1+2*cos2deltaM*cos2deltaM)-B/6*cos2deltaM*(-3+4*sin_delta*sin_delta)*(-3+4*cos2deltaM*cos2deltaM))));
    dleta_ = vpa(delta);
    delta = vpa(s / (b*A) + delta_delta);
end

tmp = vpa(sinU1*sin_delta - cosU1*cos_delta*cos_alpha1);
lat2 = vpa(atan2(sinU1*cos_delta + cosU1*sin_delta*cos_alpha1, (1-f)*sqrt(sin_alpha*sin_alpha + tmp*tmp))); %目标点纬度
lon =  vpa(atan2(sin_delta*sin_alpha1, cosU1*cos_delta - sinU1*sin_delta*cos_alpha1));
C = vpa(f/16*cosSq_alpha*(4+f*(4-3*cosSq_alpha)));
L = vpa(lon - (1-C) * f * sin_alpha *(delta + C*sin_delta*(cos2deltaM+C*cos_delta*(-1+2*cos2deltaM*cos2deltaM))));
lon2 =vpa(rem((lon1+L+3*pi) ,(2*pi)) - pi); 
% lon2 =roundn(rem((lon1+L+3*pi) ,(2*pi)) - pi, -10); %normalise to -180...+180 目标点精度

lonOut = vpa(lon2*180/pi);
latOut = vpa(lat2*180/pi);

revAz = atan2(sin_alpha, -tmp); %最终方位角

end