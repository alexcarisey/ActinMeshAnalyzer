function [ConvexAreaToArea,PerimeterAreaIndex,ContourArea,ConvexArea,Perimeter] = CheckContour(Contour)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[ConvexAreaToArea,ContourArea,ConvexArea]=ConvexAreaToContourArea(Contour);
Perimeter=ComputePerimeterLength(Contour); 
PerimeterAreaIndex=Perimeter.^2./(4.*pi*ContourArea);
Koko=size(Contour);
%if Contour(1,:)~=Contour(Koko(1),:)
%    Contour(Koko(1)+1,:)=Contour(1,:);
%end
%ImageSize=[max(max(Contour))+1 max(max(Contour))+1];
%ContourImage=CurveToImage(Contour,ImageSize);
%LabelList=bwconncomp(ContourImage);
%LabelImage=labelmatrix(LabelList);
%Stats=regionprops(LabelImage,'Eccentricity');
%

end

