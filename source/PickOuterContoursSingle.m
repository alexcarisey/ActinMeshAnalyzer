function [OuterContours,InnerContours] = PickOuterContoursSingle(Contours)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
Koko=max(size(Contours));
ContoursToRemove=[];
OuterContours=[];
InnerContours=[];
k=1;
l=1;
for i=1:Koko
    InsideSum=0;
    for j=1:Koko
            InsideSum=InsideSum+sum(inpolygon(Contours{i}(1,1),Contours{i}(1,2),Contours{j}(:,1),Contours{j}(:,2)));
    end
    InsideSum=InsideSum-1; %sum(inpolygon(Contours{i}(:,1),Contours{i}(:,2),Contours{i}(:,1),Contours{i}(:,2)));
    if InsideSum>0
        ContoursToRemove(max(size(ContoursToRemove))+1)=i;
        InnerContours{l}=Contours{i};
        l=l+1;
    else
        OuterContours{k}=Contours{i};
        k=k+1;
    end
end
NumberToRemove=max(size(ContoursToRemove));


            


end

