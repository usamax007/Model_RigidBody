classdef BodyPart < handle
    % All quantities are measured in Matlab frame
    % Each body part has an associated coordinate frame, origin and center
    % of mass which transform with the part
    % Both part and coordinate frame can be selected in the constructor 
    % to be made visible or invisible in the constructor
    % Feature of assigning children to a part is also included. Which means
    % When a part is transformed, all its children also transform with it.
    % Children must also be BodyPart objects.
    properties
        posPartOrigin; % in the Matlab frame. Not necessarily part COM position
        matPartAxes; % 3 columns representing x, y and z axes in the Matlab frame
        matColorAxes = eye(3,3);
        partName = strings
        posPartCOM;  % pos of the center of mass of the part assuming uniform density
        dataObjX;
        dataObjY;
        dataObjZ;
        objPartPlot;
        objAxesXPlot;
        objAxesYPlot;
        objAxesZPlot;
        objOriginPlot;
        objsAll; % array including handles of drawn part, drawn origin and drawn coordinate axes  
        nObjs; 
        objsChildren; % array including handles of children
        nChildren;
    end
    methods
        % constructor
        function obj = BodyPart(objType,name,dataX,dataY,dataZ,posOrigin,matAxes,children,originMarkerSize,colorValues,partVisible,axesVisible)
            obj.partName = name;
            obj.posPartOrigin = posOrigin;
            obj.matPartAxes = matAxes;
            obj.objsChildren = children;
            obj.nChildren = length(children);
            % objType is 'S' or 'V'. V for closed surface (currently only ellipsoid) and S for open surface.
            if objType == 'V' % part similar to body (head+thorax+abdomen combined) is a closed surface and modeled as an ellipsoid
                [obj.dataObjX,obj.dataObjY,obj.dataObjZ] = ellipsoid(obj.posPartOrigin(1),obj.posPartOrigin(2),obj.posPartOrigin(3),...
                    dataX,dataY,dataZ);
                obj.objPartPlot = surf(obj.dataObjX,obj.dataObjY,obj.dataObjZ,'FaceColor',colorValues,'EdgeColor', 'None',Visible=partVisible);
                % calculate part center of mass
                obj.posPartCOM = [mean(mean(obj.objPartPlot.XData));mean(mean(obj.objPartPlot.YData));mean(mean(obj.objPartPlot.ZData))];
            elseif objType == 'S' % part similar to a wing (flat plate) must be a surface and drawn as a patch or a polygon but it has a mass
                obj.dataObjX = obj.posPartOrigin(1) + dataX;
                obj.dataObjY = obj.posPartOrigin(2) + dataY;
                obj.dataObjZ = obj.posPartOrigin(3) + dataZ;
                obj.objPartPlot = patch('XData',obj.dataObjX,'YData',obj.dataObjY,'ZData',obj.dataObjZ,'FaceColor',colorValues,'EdgeColor','r',Visible=partVisible);
                obj.posPartCOM = [mean(mean(obj.objPartPlot.XData));mean(mean(obj.objPartPlot.YData));0];
            else
                error('Invalid objType');
            end
            % plot origin as a small sphere
            [originDataX,originDataY,originDataZ] = ellipsoid(obj.posPartOrigin(1),obj.posPartOrigin(2),obj.posPartOrigin(3),...
                originMarkerSize,originMarkerSize,originMarkerSize);
            obj.objOriginPlot = surf(originDataX,originDataY,originDataZ,'FaceColor','g','EdgeColor', 'None',Visible=axesVisible);
            
            % plot x-y-z axes as three orthogonal lines
            obj.objAxesXPlot = plot3(obj.posPartOrigin(1)+[0 obj.matPartAxes(1,1)],obj.posPartOrigin(2)+[0 obj.matPartAxes(2,1)],obj.posPartOrigin(3)+[0 obj.matPartAxes(3,1)],'color',obj.matColorAxes(1,:),Visible=axesVisible);
            obj.objAxesYPlot = plot3(obj.posPartOrigin(1)+[0 obj.matPartAxes(1,2)],obj.posPartOrigin(2)+[0 obj.matPartAxes(2,2)],obj.posPartOrigin(3)+[0 obj.matPartAxes(3,2)],'color',obj.matColorAxes(2,:),Visible=axesVisible);
            obj.objAxesZPlot = plot3(obj.posPartOrigin(1)+[0 obj.matPartAxes(1,3)],obj.posPartOrigin(2)+[0 obj.matPartAxes(2,3)],obj.posPartOrigin(3)+[0 obj.matPartAxes(3,3)],'color',obj.matColorAxes(3,:),Visible=axesVisible);
            
            % include all associated objects (part, origin sphere, axes lines)
            % into an array so everything is transformed together when the BodyPart is transformed 
            obj.objsAll = [obj.objPartPlot obj.objOriginPlot obj.objAxesXPlot obj.objAxesYPlot obj.objAxesZPlot]; % indices > 2 has axis which is not translated only rotated
            obj.nObjs = length(obj.objsAll);
        end
        function obj = partTransform(obj,axisRotat,angleRotat,posOriginRotat,vecTrans) % all in Matlab frame 
            % Give vectors column vectors and angle in radians
            % axisRotat is the 3x1 axis of rotation
            % angleRotat is the angle of rotation in radians
            % posOriginRotat is the position where the axis of rotation is located
            % vecTrans is the 3x1 linear translation vector 
            
            % first transform all children through recurrence
            for i_child = 1:obj.nChildren
                obj.objsChildren(i_child) = obj.objsChildren(i_child).partTransform(axisRotat,angleRotat,posOriginRotat,vecTrans);
            end

            % first rotate the entire object
            rotate(obj.objsAll,axisRotat',180/pi*angleRotat,posOriginRotat);
            % then translate the entire object
            for i_obj = 1:obj.nObjs
                    obj.objsAll(i_obj).XData = obj.objsAll(i_obj).XData + vecTrans(1);
                    obj.objsAll(i_obj).YData = obj.objsAll(i_obj).YData + vecTrans(2);
                    obj.objsAll(i_obj).ZData = obj.objsAll(i_obj).ZData + vecTrans(3);
            end
            % update origin and center of mass position
            obj.posPartOrigin = AxAng2RotMat(axisRotat,angleRotat)*(obj.posPartOrigin-posOriginRotat) + posOriginRotat + vecTrans; 
            obj.posPartCOM = AxAng2RotMat(axisRotat,angleRotat)*(obj.posPartCOM-posOriginRotat) + posOriginRotat + vecTrans;

            % update coordinate frame axes
            obj.matPartAxes = [diff(obj.objAxesXPlot.XData) diff(obj.objAxesYPlot.XData) diff(obj.objAxesZPlot.XData);...
                               diff(obj.objAxesXPlot.YData) diff(obj.objAxesYPlot.YData) diff(obj.objAxesZPlot.YData);...
                               diff(obj.objAxesXPlot.ZData) diff(obj.objAxesYPlot.ZData) diff(obj.objAxesZPlot.ZData)];
            obj.dataObjX = obj.objPartPlot.XData;
            obj.dataObjY = obj.objPartPlot.YData;
            obj.dataObjZ = obj.objPartPlot.ZData; 
        end
        function obj = updateChildren(obj,newChildren)
            obj.objsChildren = newChildren;
            obj.nChildren = length(obj.objsChildren);
        end
    end
end