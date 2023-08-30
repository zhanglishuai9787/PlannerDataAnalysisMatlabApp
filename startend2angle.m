function [yaw,s] = startend2angle(lon1,lat1,lon2,lat2)
%根据p1和p2经纬度估算两点的距离和方位角(p1 → p2)
%输入：lon1,lat1,lon2,lat2 两点经纬度
%输出：yaw 两点连线方位角deg 沿1——2方向，s两点距离 m

% lon1 = deg2rad(116.5123724) ; lat1 = deg2rad(39.7914833);  %p1	
% lon2 = deg2rad(116.5123985); lat2 = deg2rad(39.7914958);   %p2
lon1 = deg2rad(lon1) ; lat1 = deg2rad(lat1);  %p1	
lon2 = deg2rad(lon2); lat2 = deg2rad(lat2);   %p2
f = 1 /  298.257223563;
a= 6378137.0;	
b= 6356752.314245;

L = lon2 - lon1;
tanU1 = (1-f)*tan(lat1); cosU1 = 1 / sqrt((1 + tanU1*tanU1));sinU1 = tanU1 * cosU1;
tanU2 = (1-f)*tan(lat2); cosU2 = 1 / sqrt((1 + tanU2*tanU2));sinU2 = tanU2 * cosU2;
lambda = L;
lambda_ = 0;
iterationLimit = 100;
while (abs(lambda - lambda_) > 1e-12 && iterationLimit>0)
        iterationLimit = iterationLimit -1;
        sinlambda = sin(lambda);
        coslambda = cos(lambda);
        sinSq_delta = (cosU2*sinlambda) * (cosU2*sinlambda) + (cosU1*sinU2-sinU1*cosU2* coslambda) * (cosU1*sinU2-sinU1*cosU2* coslambda);
        sin_delta = sqrt(sinSq_delta);
        if sin_delta==0 
               return 
        end
        cos_delta = sinU1*sinU2 + cosU1*cosU2*coslambda;
        delta = atan2(sin_delta, cos_delta);
        sin_alpha = cosU1 * cosU2 * sinlambda / sin_delta;
        cosSq_alpha = 1 - sin_alpha*sin_alpha;
        cos2_deltaM = cos_delta - 2*sinU1*sinU2/cosSq_alpha;
        C = f/16*cosSq_alpha*(4+f*(4-3*cosSq_alpha));
        lambda_ = lambda;
        lambda = L + (1-C) * f * sin_alpha * (delta + C*sin_delta*(cos2_deltaM+C*cos_delta*(-1+2*cos2_deltaM*cos2_deltaM)));
end
uSq = cosSq_alpha * (a*a - b*b) / (b*b);
A = 1 + uSq/16384*(4096+uSq*(-768+uSq*(320-175*uSq)));
B = uSq/1024 * (256+uSq*(-128+uSq*(74-47*uSq)));
delta_delta = B*sin_delta*(cos2_deltaM+B/4*(cos_delta*(-1+2*cos2_deltaM*cos2_deltaM)-B/6*cos2_deltaM*(-3+4*sin_delta*sin_delta)*(-3+4*cos2_deltaM*cos2_deltaM)));
s = b*A*(delta-delta_delta);
fwdAz = atan2(cosU2*sinlambda,  cosU1*sinU2-sinU1*cosU2*coslambda); %初始方位角
revAz = atan2(cosU1*sinlambda, -sinU1*cosU2+cosU1*sinU2*coslambda); %最终方位角
yaw = rad2deg(fwdAz);

end