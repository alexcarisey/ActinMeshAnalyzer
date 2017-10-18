function [AreaRatio,ContourArea,ConvexArea] = ConvexAreaToContourArea(Contour)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    ContourArea=polyarea(Contour(:,1),Contour(:,2));
    DelCont=DelaunayTri(Contour(:,1),Contour(:,2));
	[~,ConvexArea]=convexHull(DelCont);
    AreaRatio=ConvexArea./ContourArea;

end

