clc
clear
close all
tic
boundary_table = readgeotable('yzsfq_hd_boundary_v22_84.shp');
crossinfo_table = readgeotable('yzsfq_hd_crossinfo_v22_84.shp');
intersection_polygon_table = readgeotable('yzsfq_hd_intersection_polygon_v22_84.shp');
lane_markingline_table = readgeotable('yzsfq_hd_lane_markingline_v22_84.shp');
lane_nodepoint_table = readgeotable('yzsfq_hd_lane_nodepoint_v22_84.shp');
laneline_table = readgeotable('yzsfq_hd_laneline_v22_84.shp');
link_nodepoint_table = readgeotable('yzsfq_hd_link_nodepoint_v22_84.shp');
linkline_table = readgeotable('yzsfq_hd_linkline_v22_84.shp');
stop_line_table = readgeotable('yzsfq_hd_stop_line_v22_84.shp');
crosswalk_table = readgeotable('yzsfq_hd_crosswalk_v22_84.shp');
save('TableOfYizhuang.mat',"boundary_table","crossinfo_table","intersection_polygon_table","lane_markingline_table","lane_nodepoint_table",...
    "laneline_table","link_nodepoint_table","linkline_table","stop_line_table","crosswalk_table")
toc
%  Maplane = shaperead('laneline_v22_84.shp');
tic
figure(1)
g_laneline =geoplot(laneline_table,Color="#D95319");
g_laneline =geoplot(laneline_table);


hold on

g_lane_nodepoint = geoplot(lane_nodepoint_table,MarkerSize=20);
g_lane_markingline =geoplot(lane_markingline_table,Color=[1,0.8,0.2]);

g_stop_line =geoplot(stop_line_table,Color=[1,0,0]);
g_intersection_polygon = geoplot(intersection_polygon_table(1,:),FaceColor=[0.8500 0.3250 0.0980]);
g_crosswalk = geoplot(crosswalk_table,FaceColor=[0.6350 0.0780 0.1840]);
g_crossinfo = geoplot(crossinfo_table,MarkerSize=20);
t_crossinfo = text(g_crossinfo.ShapeData.Latitude, g_crossinfo.ShapeData.Longitude,cellstr(crossinfo_table.cross_code)');


g_linkline = geoplot(linkline_table,Color=[0.8500 0.3250 0.0980]);
g_link_nodepoint = geoplot(link_nodepoint_table,MarkerSize=20);
g_boundary = geoplot(boundary_table,Color=[1,0.5,0.2]);
toc


lat = [39 45 19 39];
lon = [-113 -49 -100 -113];
polyshp = geopolyshape(lat,lon);

geoplot(polyshp)
hold on
tic
lat = [38 45 19 38];
lon = [-113 -49 -100 -113];
polyshp = geopolyshape(lat,lon);

h = geoplot(polyshp)
toc

datetime(1687745720264/1000, 'ConvertFrom', 'posixtime' ,'TimeZone','local','Format','uuuu-MM-dd''T''HH:mm:ss.SSS')