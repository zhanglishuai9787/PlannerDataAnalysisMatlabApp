classdef DataAnalyse < handle
    properties
        path
        T
        dataIndex
        dataIndexpre
        M
        D
        DataCursorPublisher
        dcm
        D_cell
        GeoTable  %table 地图数据
        vehData_T
        vehData_runIndex
        vehData_targets
    end
    properties(Hidden)
        Map_lane
        Map_link
        h
        t
        lat
        lon
        numOfMeasureLine
    end
    methods
        function obj = DataAnalyse(filePath)
            if exist('MapOfYizhuang_Shp.mat','file') == 2
                load('MapOfYizhuang_Shp.mat','Maplane','Maplink');
            else
                curtFilePath = pwd;
                Maplane = shaperead([curtFilePath,'\mapshpfile\yizhuang_84\','laneline_v22_84.shp']) ;
                Maplink = shaperead([curtFilePath,'\mapshpfile\yizhuang_84\','linkline.shp']) ;
                save('MapOfYizhuang_Shp.mat',"Maplink","Maplane");
            end
            if exist('TableOfYizhuang.mat','file') == 2
                load('TableOfYizhuang.mat',"boundary_table","crossinfo_table","intersection_polygon_table","lane_markingline_table","lane_nodepoint_table",...
                    "laneline_table","link_nodepoint_table","linkline_table","stop_line_table","crosswalk_table","lanes_new_table");
            else
                curtFilePath = pwd;
                boundary_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_boundary_v22_84.shp']);
                crossinfo_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_crossinfo_v22_84.shp']);
                intersection_polygon_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_intersection_polygon_v22_84.shp']);
                lane_markingline_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_lane_markingline_v22_84.shp']);
                lane_nodepoint_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_lane_nodepoint_v22_84.shp']);
                laneline_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_laneline_v22_84.shp']);
                link_nodepoint_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_link_nodepoint_v22_84.shp']);
                linkline_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_linkline_v22_84.shp']);
                stop_line_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_stop_line_v22_84.shp']);
                crosswalk_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_crosswalk_v22_84.shp']);
                lanes_new_table = readgeotable([curtFilePath,'\mapshpfile\yz\','lanes.shp']);
                save('TableOfYizhuang.mat',"boundary_table","crossinfo_table","intersection_polygon_table","lane_markingline_table","lane_nodepoint_table",...
                    "laneline_table","link_nodepoint_table","linkline_table","stop_line_table","crosswalk_table","lanes_new_table")
            end
            if isempty(filePath)
                [file, path]= uigetfile({'*.csv';'*.xlsx'});
                obj.path = [path file];
            else
                file=1;
                obj.path = filePath;
            end
            if file~=0
                [~,~,expanded_name]=fileparts(obj.path);
                if strcmp(expanded_name,'.csv')
                    opts = detectImportOptions(obj.path);
                    opts = setvaropts(opts,15,"Type","char");
                    obj.T = readtable(obj.path,opts);%读取
                elseif strcmp(expanded_name,'.xlsx')
                    opts = detectImportOptions(obj.path,"Sheet",'traj_log');
                    opts = setvaropts(opts,15,"Type","char");
                    obj.T = readtable(obj.path,opts,"Sheet",'traj_log');%读取
                    obj.vehData_T = readtable(obj.path,'Sheet','run');
                end
            end

            obj.GeoTable.boundary_table = boundary_table;
            obj.GeoTable.crossinfo_table = crossinfo_table;
            obj.GeoTable.intersection_polygon_table = intersection_polygon_table;
            obj.GeoTable.lane_markingline_table = lane_markingline_table;
            obj.GeoTable.lane_nodepoint_table = lane_nodepoint_table;
            obj.GeoTable.laneline_table = laneline_table;
            obj.GeoTable.link_nodepoint_table = link_nodepoint_table;
            obj.GeoTable.linkline_table = linkline_table;
            obj.GeoTable.stop_line_table = stop_line_table;
            obj.GeoTable.crosswalk_table = crosswalk_table;
            obj.GeoTable.lanes_new_table = lanes_new_table;
            obj.Map_lane = Maplane;
            obj.Map_link = Maplink;
            obj.M = containers.Map('KeyType','char','ValueType','char');
            obj.D = dictionary;
            obj.DataCursorPublisher = CursorPublisher();
            obj.numOfMeasureLine = 0;
            obj.dataIndexpre = 0;
            obj.dataIndex = 0;
            load('VariableNames.mat','variableNames');
            for i=1:length(variableNames)
                obj.D(variableNames{1,i}) = strcat(cell2mat(variableNames(i))," = ","Not found");
            end

%             fprintf('Map of Yizhuang\n')
        end

        function geoshowdata(obj)
            numOfpoint = nnz(obj.T.GNSS_LONG);
            obj.D_cell = cell(numOfpoint,1);
            obj.lon = obj.T.GNSS_LONG(1:numOfpoint); %经度
            obj.lat = obj.T.GNSS_LAT(1:numOfpoint); %纬度
            longitudeGoal = obj.T.longitudeGoal(1);%终点
            latitudeGoal = obj.T.latitudeGoal(1);%终点
            obj.h.p=geoplot(obj.lat,obj.lon,'Color',"#4DBEEE","Marker",'o','LineStyle','none');%在地理坐标中绘制轨迹坐标
            hold on
            if longitudeGoal ~= 0
                obj.h.goal = geoplot(latitudeGoal,longitudeGoal,Marker="pentagram");
                obj.t.goal = text(latitudeGoal, longitudeGoal,'end');
            end
            geobasemap streets-light
            obj.dcm = datacursormode;
            obj.dcm.Enable = 'on';
            obj.dcm.DisplayStyle = 'window';
            obj.dcm.UpdateFcn = @displayinfo;
%             lh = addlistener(obj.DataCursorPublisher,'dataIndeChanged',@RespondToEventCursor);
            function txt = displayinfo(~,info)
                CursorInfo=getCursorInfo(obj.dcm);
                obj.dataIndexpre = obj.dataIndex;
                obj.dataIndex = CursorInfo.DataIndex;
                if isempty(obj.D_cell{obj.dataIndex,1})
                    VariableNames = obj.T.Properties.VariableNames;%元胞数组
                    for i=1:length(VariableNames)
%                     var = table2array(obj.T(obj.dataIndex,i));
                        var = obj.T{obj.dataIndex,i};
                        if iscell(var)
                            str = strcat(cell2mat(VariableNames(i))," = ",string(cell2mat(var)));
                        elseif ~isa(var,"datetime") && isnan(var)
                            str = strcat(cell2mat(VariableNames(i))," = ","NaN");
                        else
                            str = strcat(cell2mat(VariableNames(i))," = ",string(var));

                        end
%                         disp(str);
                        obj.D(VariableNames{1,i}) = str;
                        obj.D_cell{obj.dataIndex,1} = obj.D;
                    end
                else
                    obj.D = obj.D_cell{obj.dataIndex,1};
                end

%                 for i=1:length(VariableNames)
%                       var = obj.T{obj.dataIndex,i};
%                     if iscell(var)
%                         str = strcat(cell2mat(VariableNames(i))," = ",string(cell2mat(var)));
%                     elseif ~isa(var,"datetime") && isnan(var)
%                         str = strcat(cell2mat(VariableNames(i))," = ","NaN");
%                     else
%                         str = strcat(cell2mat(VariableNames(i))," = ",string(var));
% 
%                     end
%                     obj.M(VariableNames{1,i}) = str;
%                 end
                x = info.Position(1);
                y = info.Position(2);
                txt = [
                    string(['lat',' = ',num2str(x)]);
                    string(['lon',' = ',num2str(y)]);
                    ];
                 obj.DataCursorPublisher.setdataindex(obj.dataIndex);%设置当前游标索引通报事件  （M计算完成后通报）
            end
        end

        function replayer(obj,replaymode)            
            numOfpoint = nnz(obj.T.GNSS_LONG);
            CursorInfo=getCursorInfo(obj.dcm);
            if replaymode==1 && ~isempty(CursorInfo)
                k=CursorInfo.DataIndex;   
            else
                k=1;
            end
            firstpoint = k;
            while k<=numOfpoint
                obj.dataIndexpre = obj.dataIndex;
                obj.dataIndex = k;
                if isempty(obj.D_cell{obj.dataIndex,1})
                    VariableNames = obj.T.Properties.VariableNames;%元胞数组
                    for i=1:length(VariableNames)
                        var = obj.T{obj.dataIndex,i};
                        if iscell(var)
                            str = strcat(cell2mat(VariableNames(i))," = ",string(cell2mat(var)));
                        elseif ~isa(var,"datetime") && isnan(var)
                            str = strcat(cell2mat(VariableNames(i))," = ","NaN");
                        else
                            str = strcat(cell2mat(VariableNames(i))," = ",string(var));
                        end
                        obj.D(VariableNames{1,i}) = str;
                        obj.D_cell{obj.dataIndex,1} = obj.D;
                    end
                else
                    obj.D = obj.D_cell{obj.dataIndex,1};
                end
                 obj.DataCursorPublisher.setdataindex(obj.dataIndex);%设置当前游标索引通报事件  （M计算完成后通报）
                 if k==firstpoint
                 dt = datatip(obj.h.p,'DataIndex',k);
                 elseif isvalid (dt)
                     dt.DataIndex = k;
                 else
                     break
                 end
                 pause(0.2)
                 k=k+1;
            end
        end

        function mapchange(obj,mapid)
            if mapid == 0
                if exist('MapOfYizhuang_Shp.mat','file') == 2
                    load('MapOfYizhuang_Shp.mat','Maplane','Maplink');
                else
                    curtFilePath = pwd;
                    Maplane = shaperead([curtFilePath,'\mapshpfile\yizhuang_84\','laneline_v22_84.shp']) ;
                    Maplink = shaperead([curtFilePath,'\mapshpfile\yizhuang_84\','linkline.shp']) ;
                    save('MapOfYizhuang_Shp.mat',"Maplink","Maplane");
                end
                obj.Map_lane = Maplane;
                obj.Map_link = Maplink;

                if exist('TableOfYizhuang.mat','file') == 2
                    load('TableOfYizhuang.mat',"boundary_table","crossinfo_table","intersection_polygon_table","lane_markingline_table","lane_nodepoint_table",...
                        "laneline_table","link_nodepoint_table","linkline_table","stop_line_table","crosswalk_table","lanes_new_table");
                else
                    curtFilePath = pwd;
                    boundary_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_boundary_v22_84.shp']);
                    crossinfo_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_crossinfo_v22_84.shp']);
                    intersection_polygon_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_intersection_polygon_v22_84.shp']);
                    lane_markingline_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_lane_markingline_v22_84.shp']);
                    lane_nodepoint_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_lane_nodepoint_v22_84.shp']);
                    laneline_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_laneline_v22_84.shp']);
                    link_nodepoint_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_link_nodepoint_v22_84.shp']);
                    linkline_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_linkline_v22_84.shp']);
                    stop_line_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_stop_line_v22_84.shp']);
                    crosswalk_table = readgeotable([curtFilePath,'\mapshpfile\yz\','yzsfq_hd_crosswalk_v22_84.shp']);
                    lanes_new_table = readgeotable([curtFilePath,'\mapshpfile\yz\','lanes.shp']);
                    save('TableOfYizhuang.mat',"boundary_table","crossinfo_table","intersection_polygon_table","lane_markingline_table","lane_nodepoint_table",...
                        "laneline_table","link_nodepoint_table","linkline_table","stop_line_table","crosswalk_table","lanes_new_table")
                end
                obj.GeoTable.boundary_table = boundary_table;
                obj.GeoTable.crossinfo_table = crossinfo_table;
                obj.GeoTable.intersection_polygon_table = intersection_polygon_table;
                obj.GeoTable.lane_markingline_table = lane_markingline_table;
                obj.GeoTable.lane_nodepoint_table = lane_nodepoint_table;
                obj.GeoTable.laneline_table = laneline_table;
                obj.GeoTable.link_nodepoint_table = link_nodepoint_table;
                obj.GeoTable.linkline_table = linkline_table;
                obj.GeoTable.stop_line_table = stop_line_table;
                obj.GeoTable.crosswalk_table = crosswalk_table;
                obj.GeoTable.lanes_new_table = lanes_new_table;
%                 fprintf('Map of Yizhuang\n')
            elseif mapid == 1
                if exist('MapOfChongqing_Shp.mat','file') == 2
                    load('MapOfChongqing_Shp.mat','Maplane','Maplink');
                else
                    curtFilePath = pwd;
                    Maplane = shaperead([curtFilePath,'\mapshpfile\chongqing_84\','laneline_crs84.shp']) ;
                    Maplink = shaperead([curtFilePath,'\mapshpfile\chongqing_84\','linkline_crs84.shp']) ;
                    save('MapOfChongqing_Shp.mat',"Maplink","Maplane");
                end
                obj.Map_lane = Maplane;
                obj.Map_link = Maplink;

                if exist('TableOfChongqing.mat','file') == 2
                    load('TableOfChongqing.mat',"crossinfo_table","intersection_polygon_table","lane_markingline_table","lane_nodepoint_table",...
                        "laneline_table","link_nodepoint_table","linkline_table","stop_line_table","crosswalk_table");
                else
                    curtFilePath = pwd;
                    crossinfo_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_crossinfo_crs84.shp']);
                    intersection_polygon_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_intersection_polygon_crs84.shp']);
                    lane_markingline_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_lane_markingline_crs84.shp']);
                    lane_nodepoint_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_lane_nodepoint_crs84.shp']);
                    laneline_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_laneline_crs84.shp']);
                    link_nodepoint_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_link_nodepoint_crs84.shp']);
                    linkline_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_linkline_crs84.shp']);
                    stop_line_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_stop_line_crs84.shp']);
                    crosswalk_table = readgeotable([curtFilePath,'\mapshpfile\cq\','yzsfq_hd_crosswalk_crs84.shp']);
                    save('TableOfChongqing.mat',"crossinfo_table","intersection_polygon_table","lane_markingline_table","lane_nodepoint_table",...
                        "laneline_table","link_nodepoint_table","linkline_table","stop_line_table","crosswalk_table")
                end

                obj.GeoTable.crossinfo_table = crossinfo_table;
                obj.GeoTable.intersection_polygon_table = intersection_polygon_table;
                obj.GeoTable.lane_markingline_table = lane_markingline_table;
                obj.GeoTable.lane_nodepoint_table = lane_nodepoint_table;
                obj.GeoTable.laneline_table = laneline_table;
                obj.GeoTable.link_nodepoint_table = link_nodepoint_table;
                obj.GeoTable.linkline_table = linkline_table;
                obj.GeoTable.stop_line_table = stop_line_table;
                obj.GeoTable.crosswalk_table = crosswalk_table;
%                 fprintf('Map of Chongqing\n')

            end

        end
        
        function lanelist(obj)
            if ~isempty(obj.dataIndex)
                if sum(strcmp('lanelist',fieldnames(obj.h)))~=0%存在lanelist
                    delete(obj.h.lanelist)% 删除上次画图
                    delete(obj.t.lanelist)% 删除上次画图
                end
%                 laneList = str2num(cell2mat(obj.T.laneList(obj.dataIndex)));
                laneList = jsondecode(cell2mat(obj.T.laneList(obj.dataIndex)));
                numoflanelist=length(laneList);
                index=1;
                mapOfLaneList(1:numoflanelist,1) = obj.Map_lane(1);
                for i=1:length(obj.Map_lane)
                    if sum(str2double(obj.Map_lane(i).lane_id) == laneList)
                        mapOfLaneList(index) = obj.Map_lane(i);
                        if index == numoflanelist
                            break
                        end
                        index=index+1;
                    end
                end
                % geoshow(mapOfLaneList)
                for j=1:numoflanelist
                    lonoflane=extractfield(mapOfLaneList(j),'X')';
                    latoflane=extractfield(mapOfLaneList(j),'Y')';
                    lonoflane=lonoflane(~isnan(lonoflane));
                    latoflane=latoflane(~isnan(latoflane));
                    obj.h.lanelist(j) = geoplot(latoflane,lonoflane,'LineWidth',2);
                    str = mapOfLaneList(j).lane_id;
                    obj.t.lanelist(j) = text(latoflane(1),lonoflane(1),str);
                end
            end
        end

        function linklist(obj)
            if ~isempty(obj.dataIndex)
                if sum(strcmp('linklist',fieldnames(obj.h)))~=0
                    delete(obj.h.linklist)% 删除上次画图
                    delete(obj.t.linklist)
                end
                %     Map_link = shaperead('linkline.shp') ;
                % geoshow(Map_link);%把地图画出来
                % linkList = cell2mat(T.linkList(1));
                linkList = cell2mat(obj.T.linkList(obj.dataIndex));
                linkList = strrep(linkList,char(39),'');
%                 linkList = str2num(linkList);
                linkList = jsondecode(linkList);
                numoflanelist=length(linkList);
                index=1;
                mapOfLaneList(1:numoflanelist,1) = obj.Map_link(1);
                for i=1:length(obj.Map_link)
                    if sum(str2double(obj.Map_link(i).link_id) == linkList)
                        mapOfLaneList(index) = obj.Map_link(i);
                        if index == numoflanelist
                            break
                        end
                        index=index+1;
                    end
                end
                % geoshow(mapOfLaneList)
                for j=1:numoflanelist
                    lonoflane=extractfield(mapOfLaneList(j),'X')';
                    latoflane=extractfield(mapOfLaneList(j),'Y')';
                    lonoflane=lonoflane(~isnan(lonoflane));
                    latoflane=latoflane(~isnan(latoflane));
                    obj.h.linklist(j) = geoplot(latoflane,lonoflane,'LineWidth',2.5);
                    % hold on
                    str = mapOfLaneList(j).link_id;
                    obj.t.linklist(j) = text(latoflane(1),lonoflane(1),str);
                end
            end
        end

        function linkset(obj)
            if ~isempty(obj.dataIndex)
                hold on
                if sum(strcmp('linkset',fieldnames(obj.h)))~=0
                    delete(obj.h.linkset)% 删除上次画图
                    delete(obj.t.linkset)
                end
                %     Map_link = shaperead('linkline.shp') ;
                % geoshow(Map_link);%把地图画出来
                linkSet = cell2mat(obj.T.linkSet(obj.dataIndex));
                linkSet = strrep(linkSet,char(39),'');
                linkSet = jsondecode(linkSet);
                numoflanelist=length(linkSet);
                index=1;
                mapOfLaneList(1:numoflanelist,1) = obj.Map_link(1);
                for i=1:length(obj.Map_link)
                    if sum(str2double(obj.Map_link(i).link_id) == linkSet)
                        mapOfLaneList(index) = obj.Map_link(i);
                        if index == numoflanelist
                            break
                        end
                        index=index+1;
                    end
                end
                % geoshow(mapOfLaneList)
                for j=1:numoflanelist
                    lonoflane=extractfield(mapOfLaneList(j),'X')';
                    latoflane=extractfield(mapOfLaneList(j),'Y')';
                    lonoflane=lonoflane(~isnan(lonoflane));
                    latoflane=latoflane(~isnan(latoflane));
                    obj.h.linkset(j) = geoplot(latoflane,lonoflane,'LineWidth',4);
                    % hold on
                    str = mapOfLaneList(j).link_id;
                    obj.t.linkset(j) = text(latoflane(1),lonoflane(1),str);
                end
            end
        end
        
        function curtpath(obj)
            if ~isempty(obj.dataIndex)
                if sum(strcmp('curtPath2Points',fieldnames(obj.h)))~=0%存在lanelist
                    delete(obj.h.curtPath2Points)% 删除上次画图
                end
                for index=obj.dataIndex:-1:1
                    curtPath2Points = jsondecode(cell2mat(obj.T.curtPath2Points(index)));
                    if ~isempty(curtPath2Points)
                        break
                    end
                end
                obj.h.curtPath2Points = geoplot(curtPath2Points(:,2),curtPath2Points(:,1),'r',"Marker","x",'LineStyle','-','LineWidth',2.5);
            end
        end

        function tagtPath(obj)
            if ~isempty(obj.dataIndex)
                if sum(strcmp('tagtPath2Points',fieldnames(obj.h)))~=0%存在lanelist
                    delete(obj.h.tagtPath2Points)% 删除上次画图
                end
                for index=obj.dataIndex:-1:1
                    tagtPath2Points = jsondecode(cell2mat(obj.T.tagtPath2Points(index)));
                    if ~isempty(tagtPath2Points)
                        obj.h.tagtPath2Points = geoplot(tagtPath2Points(:,2),tagtPath2Points(:,1),'r',"Marker","*",'LineStyle','-','LineWidth',2.5);
                        break
                    end
                end
            end
        end

        function envirvehinfo(obj)
            envirvehinfo = zeros(4,2);
            text_str = ["LF";"RF";"LR";"RR"];
            if ~isempty(obj.dataIndex)
                envirvehinfo(1,1) = obj.T.leftFrontLat(obj.dataIndex);
                envirvehinfo(1,2) = obj.T.leftFrontLon(obj.dataIndex);
                envirvehinfo(2,1) = obj.T.rightFrontLat(obj.dataIndex);
                envirvehinfo(2,2) = obj.T.rightFrontLon(obj.dataIndex);
                envirvehinfo(3,1) = obj.T.leftBehindLat(obj.dataIndex);
                envirvehinfo(3,2) = obj.T.leftBehindLon(obj.dataIndex);
                envirvehinfo(4,1) = obj.T.rightBehindLat(obj.dataIndex);
                envirvehinfo(4,2) = obj.T.rightBehindLon(obj.dataIndex);
                for index = 1:4
                    if envirvehinfo(index,1) ~= 0
                        obj.h.envirvehinfo(index) = geoplot(envirvehinfo(index,1),envirvehinfo(index,2),'r',"Marker","*",'LineStyle','none');
                        obj.t.envirvehinfo(index) = text(envirvehinfo(index,1),envirvehinfo(index,2),text_str(index));
                    end
                end
            end
        end

        function envvehfromavomain(obj)
            if sum(strcmp('envvehfromavomain',fieldnames(obj.h)))~=0%存在lanelist
                delete(obj.h.envvehfromavomain);% 删除上次画图
                delete(obj.t.envvehfromavomain);
            end
            envirvehinfo = zeros(8,2);
            text_str = ["F1";"F2";"F3";"F4";"R1";"R2";"R3";"R4"];
            if ~isempty(obj.dataIndex) && ~isempty(obj.T.AvoMainRoVehInfo_targetLaneFrontPos{obj.dataIndex,1})
                envirvehinfo(1:4,:) = jsondecode(obj.T.AvoMainRoVehInfo_targetLaneFrontPos{obj.dataIndex,1});
                envirvehinfo(5:8,:) = jsondecode(obj.T.AvoMainRoVehInfo_targetLaneBehindPos{obj.dataIndex,1});
%                 envirvehinfo(1,1) = obj.T.leftFrontLat(obj.dataIndex);
%                 envirvehinfo(1,2) = obj.T.leftFrontLon(obj.dataIndex);
%                 envirvehinfo(2,1) = obj.T.rightFrontLat(obj.dataIndex);
%                 envirvehinfo(2,2) = obj.T.rightFrontLon(obj.dataIndex);
%                 envirvehinfo(3,1) = obj.T.leftBehindLat(obj.dataIndex);
%                 envirvehinfo(3,2) = obj.T.leftBehindLon(obj.dataIndex);
%                 envirvehinfo(4,1) = obj.T.rightBehindLat(obj.dataIndex);
%                 envirvehinfo(4,2) = obj.T.rightBehindLon(obj.dataIndex);
                for index = 1:8
                    if envirvehinfo(index,1) ~= 0
                        obj.h.envvehfromavomain(index) = geoplot(envirvehinfo(index,2),envirvehinfo(index,1),'r',"Marker","square",'LineStyle','none');
                        obj.t.envvehfromavomain(index) = text(envirvehinfo(index,2),envirvehinfo(index,1),text_str(index));
                    end
                end
            end
        end

        function envvehfromoncoming(obj)
            if sum(strcmp('envvehfromoncoming',fieldnames(obj.h)))~=0%存在lanelist
                delete(obj.h.envvehfromoncoming);% 删除上次画图
                delete(obj.t.envvehfromoncoming);
            end
            envirvehinfo = zeros(12,2);
            text_str = ["L1_F";"L2_F";"L3_F";"L4_F";"L5_F";"L6_F";"L1_R";"L2_R";"L3_R";"L4_R";"L5_R";"L6_R"];
            if ~isempty(obj.dataIndex) && ~isempty(obj.T.AvoOncomingVehInfo_frontVehPos{obj.dataIndex,1})
                envirvehinfo(1:6,:) = jsondecode(obj.T.AvoOncomingVehInfo_frontVehPos{obj.dataIndex,1});
                envirvehinfo(7:12,:) = jsondecode(obj.T.AvoOncomingVehInfo_behindVehPos{obj.dataIndex,1});
                for index = 1:12
                    if envirvehinfo(index,1) ~= 0
                        obj.h.envvehfromoncoming(index) = geoplot(envirvehinfo(index,2),envirvehinfo(index,1),'r',"Marker","square",'LineStyle','none');
                        obj.t.envvehfromoncoming(index) = text(envirvehinfo(index,2),envirvehinfo(index,1),text_str(index));
                    end
                end
            end
        end

        function envpedfromavo(obj)
            if sum(strcmp('envpedfromavo',fieldnames(obj.h)))~=0%存在lanelist
                delete(obj.h.envpedfromavo);% 删除上次画图
                delete(obj.t.envpedfromavo);
            end
            envirpedinfo = zeros(40,2);
            a=1:40;
            text_str = string(a');
            if ~isempty(obj.dataIndex) && ~isempty(obj.T.pos_ped{obj.dataIndex,1})
                envirpedinfo(1:40,:) = jsondecode(obj.T.pos_ped{obj.dataIndex,1});
                for index = 1:40
                    if envirpedinfo(index,1) ~= 0
                        obj.h.envpedfromavo(index) = geoplot(envirpedinfo(index,2),envirpedinfo(index,1),'r',"Marker","diamond",'LineStyle','none');
                        obj.t.envpedfromavo(index) = text(envirpedinfo(index,2),envirpedinfo(index,1),text_str(index));
                    end
                end
            end
        end

        function envobjectfromveh(obj)
            if sum(strcmp('envobjectfromveh',fieldnames(obj.h)))~=0%存在lanelist
                delete(obj.h.envobjectfromveh);% 删除上次画图
                delete(obj.t.envobjectfromveh);
            end
            %获取上报时间的索引
            timeStampSet = obj.T.timeStampSet(obj.dataIndex);
            run_Index = find(obj.vehData_T.report_time == timeStampSet);
            obj.vehData_runIndex = run_Index;
            if isempty(run_Index)


            else
                veh_detection_cell = obj.vehData_T.detection(run_Index);
                veh_detection_struct = jsondecode(veh_detection_cell{1,1});
                obj.vehData_targets = veh_detection_struct;
                a=1:veh_detection_struct.targetsNum;
                text_str = string(a');
                if ~isempty(obj.dataIndex) && veh_detection_struct.targetsNum>0
                    for index = 1:veh_detection_struct.targetsNum
                        obj.h.envobjectfromveh(index) = geoplot(veh_detection_struct.targets(index).latitude,veh_detection_struct.targets(index).longitude,'b',"Marker","*",'LineStyle','none');
                        obj.t.envobjectfromveh(index) = text(veh_detection_struct.targets(index).latitude,veh_detection_struct.targets(index).longitude,text_str(index),'Color',[1,0,0]);
                    end
                end

            end

        end

        function trajectory(obj)
            if ~isempty(obj.dataIndex)
                if sum(strcmp('trajectory',fieldnames(obj.h)))~=0%存在lanelist
                    delete(obj.h.trajectory)% 删除上次画图
                end

                if ~isempty(cell2mat(obj.T.latitudeSet(obj.dataIndex)))
                    latitudeSet = jsondecode(cell2mat(obj.T.latitudeSet(obj.dataIndex)));
                    longitudeSet = jsondecode(cell2mat(obj.T.longitudeSet(obj.dataIndex)));
                    obj.h.trajectory = geoplot(latitudeSet,longitudeSet,'r.');
                end
            end
        end

        function opentext(obj,str)
            % VariableName = 'VEH_SPD';
            VariableName = str;
            delete(findobj('type','text'));
            for j=1:length(obj.lat)
%                 value_str = num2str(table2array(obj.T(j,VariableName)));
                var = table2array(obj.T(j,VariableName));
                if iscell(var)
                    value_str = string(cell2mat(var));
                elseif ~isa(var,"datetime") && isnan(var)
                    value_str = "NaN";
                else
                    value_str = string(var);
                end

                obj.t.text(j) = text(obj.lat(j),obj.lon(j),value_str);
            end
            uistack(obj.h.p,'top');
        end

        function deletealltext(obj)
            if ~isempty(obj.h)
                h_fieldnames = fieldnames(obj.h);
                if sum(strcmp('lanelist',h_fieldnames)) ~= 0
                    delete(obj.h.lanelist)
                end
                if sum(strcmp('linkset',h_fieldnames)) ~= 0
                    delete(obj.h.linkset)
                end
                if sum(strcmp('linklist',h_fieldnames)) ~= 0
                    delete(obj.h.linklist)
                end
                if sum(strcmp('curtPath2Points',h_fieldnames)) ~= 0
                    delete(obj.h.curtPath2Points)
                end
                if sum(strcmp('tagtPath2Points',h_fieldnames)) ~= 0
                    delete(obj.h.tagtPath2Points)
                end
                if sum(strcmp('envirvehinfo',h_fieldnames)) ~= 0
                    delete(obj.h.envirvehinfo)
                end
                if sum(strcmp('trajectory',h_fieldnames)) ~= 0
                    delete(obj.h.trajectory)
                end
                if sum(strcmp('envvehfromavomain',h_fieldnames)) ~= 0
                    delete(obj.h.envvehfromavomain)
                end
                if sum(strcmp('envvehfromoncoming',h_fieldnames)) ~= 0
                    delete(obj.h.envvehfromoncoming)
                end
                if sum(strcmp('envpedfromavo',h_fieldnames)) ~= 0
                    delete(obj.h.envpedfromavo)
                end
                if sum(strcmp('envobjectfromveh',h_fieldnames)) ~= 0
                    delete(obj.h.envobjectfromveh)
                end
%                 if sum(strcmp('laneOfMap',h_fieldnames)) ~= 0
%                     delete(obj.h.laneOfMap)
%                 end
%                 if sum(strcmp('allLaneOfMap',h_fieldnames)) ~= 0
%                     delete(obj.h.allLaneOfMap)
%                 end
%                 if sum(strcmp('allLinkOfMap',h_fieldnames)) ~= 0
%                     delete(obj.h.allLinkOfMap)
%                 end
%                 if sum(strcmp('linkOfMap',h_fieldnames)) ~= 0
%                     delete(obj.h.linkOfMap)
%                 end
            end
            if ~isempty(obj.t)
                t_fieldnames = fieldnames(obj.t);
                if sum(strcmp('lanelist',t_fieldnames)) ~= 0
                    delete(obj.t.lanelist)
                end
                if sum(strcmp('linklist',t_fieldnames)) ~= 0
                    delete(obj.t.linklist)
                end
                if sum(strcmp('linkset',t_fieldnames)) ~= 0
                    delete(obj.t.linkset)
                end
                if sum(strcmp('text',t_fieldnames)) ~= 0
                    delete(obj.t.text)
                end
                if sum(strcmp('envirvehinfo',t_fieldnames)) ~= 0
                    delete(obj.t.envirvehinfo)
                end
                if sum(strcmp('envvehfromavomain',t_fieldnames)) ~= 0
                    delete(obj.t.envvehfromavomain)
                end
                if sum(strcmp('envvehfromoncoming',t_fieldnames)) ~= 0
                    delete(obj.t.envvehfromoncoming)
                end
                if sum(strcmp('envpedfromavo',t_fieldnames)) ~= 0
                    delete(obj.t.envpedfromavo)
                end
                if sum(strcmp('envobjectfromveh',t_fieldnames)) ~= 0
                    delete(obj.t.envobjectfromveh)
                end
            end
        end
        
        function deletemaptext(obj)
            if ~isempty(obj.h)
                h_fieldnames = fieldnames(obj.h);
                if sum(strcmp('laneOfMap',h_fieldnames)) ~= 0
                    delete(obj.h.laneOfMap)
                end
                if sum(strcmp('allLaneOfMap',h_fieldnames)) ~= 0
                    delete(obj.h.allLaneOfMap)
                end
                if sum(strcmp('allLinkOfMap',h_fieldnames)) ~= 0
                    delete(obj.h.allLinkOfMap)
                end
                if sum(strcmp('linkOfMap',h_fieldnames)) ~= 0
                    delete(obj.h.linkOfMap)
                end
            end
        end

        function dispLaneOfMap(obj)
            %获取车辆轨迹第一个坐标点和最后一个坐标点
            latStart=obj.lat(1);
            lonStart=obj.lon(1);
            latEnd=obj.lat(end);
            lonEnd=obj.lon(end);
            %计算两点距离的一半及中点经纬度
            rOfArea = SphereDist([lonStart,latStart],[lonEnd,latEnd])*1000/2 + 50;
            latMiddle = (latStart+latEnd)/2;
            lonMiddle = (lonStart+lonEnd)/2;
            %去掉范围外的lane
            index2delete=zeros(1,length(obj.Map_lane));
            for i=1:length(obj.Map_lane)
                if SphereDist(obj.Map_lane(i).BoundingBox(1,:),[lonMiddle,latMiddle])*1000>rOfArea && SphereDist(obj.Map_lane(i).BoundingBox(2,:),[lonMiddle,latMiddle])*1000>rOfArea
                    index2delete(1,i)=1;
                end
            end
            Lane = obj.Map_lane;
            Lane(logical(index2delete)) = [];%删除范围外的lane
            %
            lineX = {};
            lineY = {};
            lineId = {};
            lineIndex=1;
            while ~isempty(Lane)
%                 if lineIndex == 414
%                     lineIndex
%                 end
                x = Lane(1).X(~isnan(Lane(1).X));
                y = Lane(1).Y(~isnan(Lane(1).Y));
                lineid = str2double(Lane(1).lane_id);
                suc_lanes = Lane(1).suc_lanes;
                pre_lanes = Lane(1).pre_lanes;
                Lane(1) = [];

                while ~isempty(suc_lanes)
                    sucFlag = 0;
                    suc_lanes_cell = strsplit(suc_lanes,':');
                    for i=1:length(suc_lanes_cell)
                        sucLaneIndex = find(strcmp({Lane.lane_id},suc_lanes_cell{1,i}) == 1);
                        if ~isempty(sucLaneIndex) && Lane(sucLaneIndex).turn_flag == 4%不优先选择掉头
                            continue
                        end
                        if ~isempty(sucLaneIndex)
                            x = [x,Lane(sucLaneIndex).X(~isnan(Lane(sucLaneIndex).X))];
                            y = [y,Lane(sucLaneIndex).Y(~isnan(Lane(sucLaneIndex).Y))];
                            lineid = [lineid,str2double(Lane(sucLaneIndex).lane_id)];
                            suc_lanes = Lane(sucLaneIndex).suc_lanes;
                            Lane(sucLaneIndex) = [];
                            sucFlag = 1;
                            break
                        end
                    end
                    if sucFlag == 0
                        suc_lanes=[];
                    end
                end
                while ~isempty(pre_lanes)
                    preFlag = 0;
                    pre_lanes_cell = strsplit(pre_lanes,':');
                    for i=1:length(pre_lanes_cell)
                        preLaneIndex = find(strcmp({Lane.lane_id},pre_lanes_cell{1,i}) == 1);
                        if ~isempty(preLaneIndex) && Lane(preLaneIndex).turn_flag == 4
                            continue
                        end
                        if ~isempty(preLaneIndex)
                            x = [Lane(preLaneIndex).X(~isnan(Lane(preLaneIndex).X)),x];
                            y = [Lane(preLaneIndex).Y(~isnan(Lane(preLaneIndex).Y)),y];
                            lineid = [str2double(Lane(preLaneIndex).lane_id),lineid];
                            pre_lanes = Lane(preLaneIndex).pre_lanes;
                            Lane(preLaneIndex) = [];
                            preFlag = 1;
                            break
                        end
                    end
                    if preFlag == 0
                        pre_lanes=[];
                    end
                end
                lineX{lineIndex,1} = x;
                lineY{lineIndex,1} = y;
                lineId{lineIndex,1} = lineid;
                lineIndex=lineIndex+1;
            end
            for i=1:length(lineX)
                obj.h.laneOfMap(i) = geoplot(lineY{i,1},lineX{i,1},'Color',"#D95319");
            end
            uistack(obj.h.p,'top');
        end

        function dispLaneOfMapManual(obj)
            while(1)
                obj.h.circleFromLaneMap = drawcircle;

                rOfArea = obj.h.circleFromLaneMap.Radius*pi/180*1000*6371;
                latMiddle = obj.h.circleFromLaneMap.Center(1) ;
                lonMiddle = obj.h.circleFromLaneMap.Center(2);

                %去掉范围外的lane
                index2delete=zeros(1,length(obj.Map_lane));
                for i=1:length(obj.Map_lane)
                    if SphereDist(obj.Map_lane(i).BoundingBox(1,:),[lonMiddle,latMiddle])*1000>rOfArea && SphereDist(obj.Map_lane(i).BoundingBox(2,:),[lonMiddle,latMiddle])*1000>rOfArea
                        index2delete(1,i)=1;
                    end
                end
                Lane = obj.Map_lane;
                Lane(logical(index2delete)) = [];%删除范围外的lane
                %
                lineX = {};
                lineY = {};
                lineId = {};
                lineIndex=1;
                while ~isempty(Lane)
                    %                 if lineIndex == 414
                    %                     lineIndex
                    %                 end
                    x = Lane(1).X(~isnan(Lane(1).X));
                    y = Lane(1).Y(~isnan(Lane(1).Y));
                    lineid = str2double(Lane(1).lane_id);
                    suc_lanes = Lane(1).suc_lanes;
                    pre_lanes = Lane(1).pre_lanes;
                    Lane(1) = [];

                    while ~isempty(suc_lanes)
                        sucFlag = 0;
                        suc_lanes_cell = strsplit(suc_lanes,':');
                        for i=1:length(suc_lanes_cell)
                            sucLaneIndex = find(strcmp({Lane.lane_id},suc_lanes_cell{1,i}) == 1);
                            if ~isempty(sucLaneIndex) && Lane(sucLaneIndex).turn_flag == 4%不优先选择掉头
                                continue
                            end
                            if ~isempty(sucLaneIndex)
                                x = [x,Lane(sucLaneIndex).X(~isnan(Lane(sucLaneIndex).X))];
                                y = [y,Lane(sucLaneIndex).Y(~isnan(Lane(sucLaneIndex).Y))];
                                lineid = [lineid,str2double(Lane(sucLaneIndex).lane_id)];
                                suc_lanes = Lane(sucLaneIndex).suc_lanes;
                                Lane(sucLaneIndex) = [];
                                sucFlag = 1;
                                break
                            end
                        end
                        if sucFlag == 0
                            suc_lanes=[];
                        end
                    end
                    while ~isempty(pre_lanes)
                        preFlag = 0;
                        pre_lanes_cell = strsplit(pre_lanes,':');
                        for i=1:length(pre_lanes_cell)
                            preLaneIndex = find(strcmp({Lane.lane_id},pre_lanes_cell{1,i}) == 1);
                            if ~isempty(preLaneIndex) && Lane(preLaneIndex).turn_flag == 4
                                continue
                            end
                            if ~isempty(preLaneIndex)
                                x = [Lane(preLaneIndex).X(~isnan(Lane(preLaneIndex).X)),x];
                                y = [Lane(preLaneIndex).Y(~isnan(Lane(preLaneIndex).Y)),y];
                                lineid = [str2double(Lane(preLaneIndex).lane_id),lineid];
                                pre_lanes = Lane(preLaneIndex).pre_lanes;
                                Lane(preLaneIndex) = [];
                                preFlag = 1;
                                break
                            end
                        end
                        if preFlag == 0
                            pre_lanes=[];
                        end
                    end
                    lineX{lineIndex,1} = x;
                    lineY{lineIndex,1} = y;
                    lineId{lineIndex,1} = lineid;
                    lineIndex=lineIndex+1;
                end
                for i=1:length(lineX)
                    obj.h.laneOfMap(i) = geoplot(lineY{i,1},lineX{i,1},'Color',"#D95319");
                end

                answer = questdlg('Do you want to choose this area?', ...
                	'Lane Map', ...
                	'Sure','Reselect','Cancel','Sure');
                if isempty(answer)
                    delete(obj.h.circleFromLaneMap)
                    delete(obj.h.laneOfMap)
                    break
                else
                    if strcmp(answer,'Sure')
                        delete(obj.h.circleFromLaneMap)
                        break
                    elseif strcmp(answer,'Reselect')
                        delete(obj.h.circleFromLaneMap)
                        delete(obj.h.laneOfMap)
                    elseif strcmp(answer,'Cancel')
                        delete(obj.h.circleFromLaneMap)
                        delete(obj.h.laneOfMap)
                        break
                    end
                end
            end
            uistack(obj.h.p,'top');
        end
          
        function dispAllLaneOfMap(obj)
            Lane = obj.Map_lane;
            lineX = {};
            lineY = {};
            lineId = {};
            lineIndex=1;
            while ~isempty(Lane)
                %                 if lineIndex == 414
                %                     lineIndex
                %                 end
                x = Lane(1).X(~isnan(Lane(1).X));
                y = Lane(1).Y(~isnan(Lane(1).Y));
                lineid = str2double(Lane(1).lane_id);
                suc_lanes = Lane(1).suc_lanes;
                pre_lanes = Lane(1).pre_lanes;
                Lane(1) = [];

                while ~isempty(suc_lanes)
                    sucFlag = 0;
                    suc_lanes_cell = strsplit(suc_lanes,':');
                    for i=1:length(suc_lanes_cell)
                        sucLaneIndex = find(strcmp({Lane.lane_id},suc_lanes_cell{1,i}) == 1);
                        if ~isempty(sucLaneIndex) && Lane(sucLaneIndex).turn_flag == 4%不优先选择掉头
                            continue
                        end
                        if ~isempty(sucLaneIndex)
                            x = [x,Lane(sucLaneIndex).X(~isnan(Lane(sucLaneIndex).X))];
                            y = [y,Lane(sucLaneIndex).Y(~isnan(Lane(sucLaneIndex).Y))];
                            lineid = [lineid,str2double(Lane(sucLaneIndex).lane_id)];
                            suc_lanes = Lane(sucLaneIndex).suc_lanes;
                            Lane(sucLaneIndex) = [];
                            sucFlag = 1;
                            break
                        end
                    end
                    if sucFlag == 0
                        suc_lanes=[];
                    end
                end
                while ~isempty(pre_lanes)
                    preFlag = 0;
                    pre_lanes_cell = strsplit(pre_lanes,':');
                    for i=1:length(pre_lanes_cell)
                        preLaneIndex = find(strcmp({Lane.lane_id},pre_lanes_cell{1,i}) == 1);
                        if ~isempty(preLaneIndex) && Lane(preLaneIndex).turn_flag == 4
                            continue
                        end
                        if ~isempty(preLaneIndex)
                            x = [Lane(preLaneIndex).X(~isnan(Lane(preLaneIndex).X)),x];
                            y = [Lane(preLaneIndex).Y(~isnan(Lane(preLaneIndex).Y)),y];
                            lineid = [str2double(Lane(preLaneIndex).lane_id),lineid];
                            pre_lanes = Lane(preLaneIndex).pre_lanes;
                            Lane(preLaneIndex) = [];
                            preFlag = 1;
                            break
                        end
                    end
                    if preFlag == 0
                        pre_lanes=[];
                    end
                end
                lineX{lineIndex,1} = x;
                lineY{lineIndex,1} = y;
                lineId{lineIndex,1} = lineid;
                lineIndex=lineIndex+1;
            end
            for i=1:length(lineX)
                obj.h.allLaneOfMap(i) = geoplot(lineY{i,1},lineX{i,1},'Color',"#D95319");
                hold on
            end
            uistack(obj.h.p,'top');
        end
        
        function dispAllLinkOfMap(obj)
            Link = obj.Map_link;
            lineX = {};
            lineY = {};
            lineId = {};
            lineIndex=1;
            while ~isempty(Link)
                x = Link(1).X(~isnan(Link(1).X));
                y = Link(1).Y(~isnan(Link(1).Y));
                lineid = str2double(Link(1).link_id);
                suc_links = Link(1).suc_links;
                pre_links = Link(1).pre_links;
                Link(1) = [];

                while ~isempty(suc_links)
                    sucFlag = 0;
                    suc_links_cell = strsplit(suc_links,':');
                    for i=1:length(suc_links_cell)
                        sucLinkIndex = find(strcmp({Link.link_id},suc_links_cell{1,i}) == 1);
                        if ~isempty(sucLinkIndex) && Link(sucLinkIndex).turn_flag == 4%不优先选择掉头
                            continue
                        end
                        if ~isempty(sucLinkIndex)
                            x = [x,Link(sucLinkIndex).X(~isnan(Link(sucLinkIndex).X))];
                            y = [y,Link(sucLinkIndex).Y(~isnan(Link(sucLinkIndex).Y))];
                            lineid = [lineid,str2double(Link(sucLinkIndex).link_id)];
                            suc_links = Link(sucLinkIndex).suc_links;
                            Link(sucLinkIndex) = [];
                            sucFlag = 1;
                            break
                        end
                    end
                    if sucFlag == 0
                        suc_links=[];
                    end
                end
                while ~isempty(pre_links)
                    preFlag = 0;
                    pre_links_cell = strsplit(pre_links,':');
                    for i=1:length(pre_links_cell)
                        preLinkIndex = find(strcmp({Link.link_id},pre_links_cell{1,i}) == 1);
                        if ~isempty(preLinkIndex) && Link(preLinkIndex).turn_flag == 4
                            continue
                        end
                        if ~isempty(preLinkIndex)
                            x = [Link(preLinkIndex).X(~isnan(Link(preLinkIndex).X)),x];
                            y = [Link(preLinkIndex).Y(~isnan(Link(preLinkIndex).Y)),y];
                            lineid = [str2double(Link(preLinkIndex).link_id),lineid];
                            pre_links = Link(preLinkIndex).pre_links;
                            Link(preLinkIndex) = [];
                            preFlag = 1;
                            break
                        end
                    end
                    if preFlag == 0
                        pre_links=[];
                    end
                end
                lineX{lineIndex,1} = x;
                lineY{lineIndex,1} = y;
                lineId{lineIndex,1} = lineid;
                lineIndex=lineIndex+1;
            end
            for i=1:length(lineX)
                obj.h.allLinkOfMap(i) = geoplot(lineY{i,1},lineX{i,1},'Color',"#EDB120");
                hold on
            end
            uistack(obj.h.p,'top');
        end
        
        function dispLinkOfMap(obj)
            %获取车辆轨迹第一个坐标点和最后一个坐标点
            latStart=obj.lat(1);
            lonStart=obj.lon(1);
            latEnd=obj.lat(end);
            lonEnd=obj.lon(end);
            %计算两点距离的一半及中点经纬度
            rOfArea = SphereDist([lonStart,latStart],[lonEnd,latEnd])*1000/2 + 50;
            latMiddle = (latStart+latEnd)/2;
            lonMiddle = (lonStart+lonEnd)/2;
            %去掉范围外的lane
            index2delete=zeros(1,length(obj.Map_link));
            for i=1:length(obj.Map_link)
                if SphereDist(obj.Map_link(i).BoundingBox(1,:),[lonMiddle,latMiddle])*1000>rOfArea && SphereDist(obj.Map_link(i).BoundingBox(2,:),[lonMiddle,latMiddle])*1000>rOfArea
                    index2delete(1,i)=1;
                end
            end


            Link = obj.Map_link;
            Link(logical(index2delete)) = [];%删除范围外的Link
            lineX = {};
            lineY = {};
            lineId = {};
            lineIndex=1;
            while ~isempty(Link)
                x = Link(1).X(~isnan(Link(1).X));
                y = Link(1).Y(~isnan(Link(1).Y));
                lineid = str2double(Link(1).link_id);
                suc_links = Link(1).suc_links;
                pre_links = Link(1).pre_links;
                Link(1) = [];

                while ~isempty(suc_links)
                    sucFlag = 0;
                    suc_links_cell = strsplit(suc_links,':');
                    for i=1:length(suc_links_cell)
                        sucLinkIndex = find(strcmp({Link.link_id},suc_links_cell{1,i}) == 1);
                        if ~isempty(sucLinkIndex) && Link(sucLinkIndex).turn_flag == 4%不优先选择掉头
                            continue
                        end
                        if ~isempty(sucLinkIndex)
                            x = [x,Link(sucLinkIndex).X(~isnan(Link(sucLinkIndex).X))];
                            y = [y,Link(sucLinkIndex).Y(~isnan(Link(sucLinkIndex).Y))];
                            lineid = [lineid,str2double(Link(sucLinkIndex).link_id)];
                            suc_links = Link(sucLinkIndex).suc_links;
                            Link(sucLinkIndex) = [];
                            sucFlag = 1;
                            break
                        end
                    end
                    if sucFlag == 0
                        suc_links=[];
                    end
                end
                while ~isempty(pre_links)
                    preFlag = 0;
                    pre_links_cell = strsplit(pre_links,':');
                    for i=1:length(pre_links_cell)
                        preLinkIndex = find(strcmp({Link.link_id},pre_links_cell{1,i}) == 1);
                        if ~isempty(preLinkIndex) && Link(preLinkIndex).turn_flag == 4
                            continue
                        end
                        if ~isempty(preLinkIndex)
                            x = [Link(preLinkIndex).X(~isnan(Link(preLinkIndex).X)),x];
                            y = [Link(preLinkIndex).Y(~isnan(Link(preLinkIndex).Y)),y];
                            lineid = [str2double(Link(preLinkIndex).link_id),lineid];
                            pre_links = Link(preLinkIndex).pre_links;
                            Link(preLinkIndex) = [];
                            preFlag = 1;
                            break
                        end
                    end
                    if preFlag == 0
                        pre_links=[];
                    end
                end
                lineX{lineIndex,1} = x;
                lineY{lineIndex,1} = y;
                lineId{lineIndex,1} = lineid;
                lineIndex=lineIndex+1;
            end
            for i=1:length(lineX)
                obj.h.linkOfMap(i) = geoplot(lineY{i,1},lineX{i,1},'Color',"#EDB120");
            end
            uistack(obj.h.p,'top');
        end

        function measureline(obj)
            obj.numOfMeasureLine = obj.numOfMeasureLine+1;
            obj.h.measureline(obj.numOfMeasureLine) = drawline('Color',[0.9290 0.6940 0.1250]);
            addlistener(obj.h.measureline(obj.numOfMeasureLine),'MovingROI',@DataAnalyse.drawlineevents);
            addlistener(obj.h.measureline(obj.numOfMeasureLine),'ROIMoved',@DataAnalyse.drawlineevents);
            lonPoint1 = obj.h.measureline(obj.numOfMeasureLine).Position(1,2);
            latPoint1 = obj.h.measureline(obj.numOfMeasureLine).Position(1,1);
            lonPoint2 = obj.h.measureline(obj.numOfMeasureLine).Position(2,2);
            latPoint2 = obj.h.measureline(obj.numOfMeasureLine).Position(2,1);
%             latMiddle = (latPoint1+latPoint2)/2;
%             lonMiddle = (lonPoint1+lonPoint2)/2;
            dist = SphereDist([lonPoint1,latPoint1],[lonPoint2,latPoint2])*1000;
            str = strcat(string(dist),"m");
            obj.h.measureline(obj.numOfMeasureLine).Label = str;
%             obj.t.measureline(obj.numOfMeasureLine) = text(latMiddle,lonMiddle,str);
        end

        function deletemeasureline(obj)
            delete(obj.h.measureline);
%             delete(obj.t.measureline);
            obj.numOfMeasureLine = 0;
        end


    end
    methods(Static)
        function drawlineevents(src,evt)
            evname = evt.EventName;
            switch(evname)
                case{'MovingROI'}
%                     disp(['ROI moving previous position: ' mat2str(evt.PreviousPosition)]);
%                     disp(['ROI moving current position: ' mat2str(evt.CurrentPosition)]);
                    lonPoint1 = evt.CurrentPosition(1,2);
                    latPoint1 = evt.CurrentPosition(1,1);
                    lonPoint2 = evt.CurrentPosition(2,2);
                    latPoint2 = evt.CurrentPosition(2,1);
                    dist = SphereDist([lonPoint1,latPoint1],[lonPoint2,latPoint2])*1000;
                    str = strcat(string(dist),"m");
                    src.Label = str;
                case{'ROIMoved'}
%                     disp(['ROI moved previous position: ' mat2str(evt.PreviousPosition)]);
%                     disp(['ROI moved current position: ' mat2str(evt.CurrentPosition)]);
                    lonPoint1 = evt.CurrentPosition(1,2);
                    latPoint1 = evt.CurrentPosition(1,1);
                    lonPoint2 = evt.CurrentPosition(2,2);
                    latPoint2 = evt.CurrentPosition(2,1);
                    dist = SphereDist([lonPoint1,latPoint1],[lonPoint2,latPoint2])*1000;
                    str = strcat(string(dist),"m");
                    src.Label = str;
            end
        end

    end
end