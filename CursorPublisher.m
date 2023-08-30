% 发布者定义事件+定义触发事件
classdef CursorPublisher < handle
    properties
        dataIndex=0;
    end
    events% 定义事件
        dataIndeChanged
    end
    methods
        % 定义触发事件
        function setdataindex(obj,newVal)
            if obj.dataIndex ~= newVal
                obj.dataIndex = newVal;
                notify(obj,'dataIndeChanged');
            end
        end
    end


end

