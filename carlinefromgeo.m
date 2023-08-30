function [FRwheel,FLwheel,RRwheel,RLwheel,body,lon2,lat2] = carlinefromgeo(lon,lat,yaw)

% lon=10;lat=10;yaw=0;

gap=0.8;
wheelWidth = 0.5;
wheelDia = 0.8;
vehL = 5;
vehW = 1.5;

[lon2,lat2]=startangle2endhigh(lon,lat,yaw,0.5*vehL);
lon_vector = lon-lon2;
lat_vector = lat-lat2;
FW_lon1 = lon+(lon2-lon)*(1-2*gap/vehL);
FW_lat1 = lat+(lat2-lat)*(1-2*gap/vehL);
ll=sqrt((lon2-lon)^2+(lat2-lat)^2);



FW_lon1r = FW_lon1+ll*(wheelWidth+vehW)/vehL*sind(yaw+90);
FW_lat1r = FW_lat1+ll*(wheelWidth+vehW)/vehL*cosd(yaw+90);
FW_lon1l = 2*FW_lon1-FW_lon1r;
FW_lat1l = 2*FW_lat1-FW_lat1r;
FW_lon2r = FW_lon1r+lon_vector*2*wheelDia/vehL;
FW_lat2r = FW_lat1r+lat_vector*2*wheelDia/vehL;
FW_lon2l = FW_lon1l+lon_vector*2*wheelDia/vehL;
FW_lat2l = FW_lat1l+lat_vector*2*wheelDia/vehL;
RW_lon2r = FW_lon1r+lon_vector*(1-2*gap/vehL)*2;
RW_lat2r = FW_lat1r+lat_vector*(1-2*gap/vehL)*2;
RW_lon1r = RW_lon2r-lon_vector*2*wheelDia/vehL;
RW_lat1r = RW_lat2r-lat_vector*2*wheelDia/vehL;
RW_lon2l = FW_lon1l+lon_vector*(1-2*gap/vehL)*2;
RW_lat2l = FW_lat1l+lat_vector*(1-2*gap/vehL)*2;
RW_lon1l = RW_lon2l-lon_vector*2*wheelDia/vehL;
RW_lat1l = RW_lat2l-lat_vector*2*wheelDia/vehL;

bd_lonfr = lon2+(FW_lon1r-FW_lon1)*vehW/(vehW+wheelWidth);
bd_latfr = lat2+(FW_lat1r-FW_lat1)*vehW/(vehW+wheelWidth);
bd_lonfl = 2*lon2-bd_lonfr;
bd_latfl = 2*lat2-bd_latfr;
bd_lonrr = bd_lonfr+lon_vector*2;
bd_latrr = bd_latfr+lat_vector*2;
bd_lonrl = bd_lonfl+lon_vector*2;
bd_latrl = bd_latfl+lat_vector*2;

FW_lon3r = bd_lonfr+lon_vector*2*(wheelDia+gap)/vehL;
FW_lat3r = bd_latfr+lat_vector*2*(wheelDia+gap)/vehL;
FW_lon4r = bd_lonfr+lon_vector*2*gap/vehL;
FW_lat4r = bd_latfr+lat_vector*2*gap/vehL;
FW_lon3l = bd_lonfl+lon_vector*2*(wheelDia+gap)/vehL;
FW_lat3l = bd_latfl+lat_vector*2*(wheelDia+gap)/vehL;
FW_lon4l = bd_lonfl+lon_vector*2*gap/vehL;
FW_lat4l = bd_latfl+lat_vector*2*gap/vehL;
RW_lon3r = bd_lonrr-lon_vector*2*gap/vehL;
RW_lat3r = bd_latrr-lat_vector*2*gap/vehL;
RW_lon4r = bd_lonrr-lon_vector*2*(wheelDia+gap)/vehL;
RW_lat4r = bd_latrr-lat_vector*2*(wheelDia+gap)/vehL;
RW_lon3l = bd_lonrl-lon_vector*2*gap/vehL;
RW_lat3l = bd_latrl-lat_vector*2*gap/vehL;
RW_lon4l = bd_lonrl-lon_vector*2*(wheelDia+gap)/vehL;
RW_lat4l = bd_latrl-lat_vector*2*(wheelDia+gap)/vehL;


FRwheel = [FW_lat4r,FW_lat1r,FW_lat2r,FW_lat3r;FW_lon4r,FW_lon1r,FW_lon2r,FW_lon3r];
FLwheel = [FW_lat4l,FW_lat1l,FW_lat2l,FW_lat3l;FW_lon4l,FW_lon1l,FW_lon2l,FW_lon3l];
RRwheel = [RW_lat4r,RW_lat1r,RW_lat2r,RW_lat3r;RW_lon4r,RW_lon1r,RW_lon2r,RW_lon3r];
RLwheel = [RW_lat4l,RW_lat1l,RW_lat2l,RW_lat3l;RW_lon4l,RW_lon1l,RW_lon2l,RW_lon3l];
body = [bd_latfr,bd_latfl,bd_latrl,bd_latrr,bd_latfr;bd_lonfr,bd_lonfl,bd_lonrl,bd_lonrr,bd_lonfr];


% body = geoplot([double(bd_latfr),double(bd_latfl),double(bd_latrl),double(bd_latrr),double(bd_latfr)],[double(bd_lonfr),double(bd_lonfl),double(bd_lonrl),double(bd_lonrr),double(bd_lonfr)],'LineWidth',2);
% hold on
% FRwheel = geoplot([double(FW_lat4r),double(FW_lat1r),double(FW_lat2r),double(FW_lat3r)],[double(FW_lon4r),double(FW_lon1r),double(FW_lon2r),double(FW_lon3r)],'LineWidth',2,'Color','k');
% FLwheel = geoplot([double(FW_lat4l),double(FW_lat1l),double(FW_lat2l),double(FW_lat3l)],[double(FW_lon4l),double(FW_lon1l),double(FW_lon2l),double(FW_lon3l)],'LineWidth',2,'Color','k');
% RRwheel = geoplot([double(RW_lat4r),double(RW_lat1r),double(RW_lat2r),double(RW_lat3r)],[double(RW_lon4r),double(RW_lon1r),double(RW_lon2r),double(RW_lon3r)],'LineWidth',2,'Color','k');
% RLwheel = geoplot([double(RW_lat4l),double(RW_lat1l),double(RW_lat2l),double(RW_lat3l)],[double(RW_lon4l),double(RW_lon1l),double(RW_lon2l),double(RW_lon3l)],'LineWidth',2,'Color','k');

% 
% geoplot(lat,lon,'')
% ax = worldmap('World');
% geoshow(ax,flip([double(bd_latfr),double(bd_latfl),double(bd_latrl),double(bd_latrr),double(bd_latfr)]),flip([double(bd_lonfr),double(bd_lonfl),double(bd_lonrl),double(bd_lonrr),double(bd_lonfr)]),'DisplayType','polygon')
% geoplot(lon,lat,'*',double(lat2),double(lon2),'^',double(lat1),double(lon1),'o',...
%     double(FW_lat1),double(FW_lon1),'x',double(FW_lat2),double(FW_lon2),'x',double(RW_lat1),double(RW_lon1),'x',double(RW_lat2),double(RW_lon2),'x'...
%     ,double(ledRlat),double(ledRlon),'x',double(ledRRlat),double(ledRRlon),'x',double(ledFlat),double(ledFlon),'x',double(ledFFlat),double(ledFFlon),'x',...
%     double(FW_lat1r),double(FW_lon1r),'x',double(FW_lat1l),double(FW_lon1l),'x', double(FW_lat2r),double(FW_lon2r),'x',double(FW_lat2l),double(FW_lon2l),'x',...
%     double(RW_lat1r),double(RW_lon1r),'x',double(RW_lat1l),double(RW_lon1l),'x', double(RW_lat2r),double(RW_lon2r),'x',double(RW_lat2l),double(RW_lon2l),'x',...
%     double(bd_latfr),double(bd_lonfr),'x',double(bd_latfl),double(bd_lonfl),'x', double(bd_latrr),double(bd_lonrr),'x',double(bd_latrl),double(bd_lonrl),'x')


end