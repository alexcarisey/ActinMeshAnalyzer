function [Perimeter] = ComputePerimeterLength(Contour)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Koko=size(Contour);
Contour2=Contour;
Contour2(Koko(1)+1,:)=Contour(1,:);
Contour2(1,:)=[];
DiffSq=(Contour-Contour2).^2;
DiffS=DiffSq(:,1)+DiffSq(:,2);
Diff=DiffS.^0.5;
Perimeter=sum(Diff);

end

