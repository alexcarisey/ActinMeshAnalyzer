function [Contours,LargeContours,SmallContours,OutsideClusterLarge,PIindexSum] = ContourOrganizerWithout(ContourMatrix,LargeLimit)
%This function converts the contour matrix into cell array of contours, and sorts the contours based on size and relation to clusters and estimated cell positions    UNTITLED2 Summary of this function goes here
%  
Koko=size(ContourMatrix);
%ContourImage=zeros(ImageSize);
%ContourImage2=zeros(ImageSize);
n=0;
p=1;
q=1;
r=1;
%MaxArea=0.25.*pi.*MaxCellD.^2;
Contours=[];
%AllLargeContours=[];
%AllSmallContours=[];
LargeContours=[];
SmallContours=[];
%OutsideClusterLarge=[];
%List=[];
i=1;
j=1;
%KokoEstCell=size(EstimatedCellPositions);
%PIindexList=zeros(KokoEstCell,1);
%PIindexSum=zeros(KokoEstCell,1);
IndexList=1:1:Koko(2);
while i<Koko(2)
    ContourPoints=ContourMatrix(2,i);
    LogicalLeft=IndexList>i;
    LogicalRight=IndexList<(i+ContourPoints+1);
    LogicalComplete=LogicalLeft.*LogicalRight;
    Apu=ContourMatrix(:,LogicalComplete==1);
    Contours{j}=transpose(Apu);
    j=j+1;
    i=i+ContourPoints+1;
end
SizeContours=max(size(Contours));
%Contours

for i=1:SizeContours
    MaxXY=max(Contours{i},[],1);
    MinXY=min(Contours{i},[],1);
    DiffX=MaxXY(1)-MinXY(1);
    DiffY=MaxXY(2)-MinXY(2);
%     %if (max(DiffX,DiffY)<(4*MaxCellD))&&(max(DiffX,DiffY)<(10*min(DiffX,DiffY)))
%     %NumberAll=inpolygon(AllCellsCoord(:,1),AllCellsCoord(:,2),Contours{i}(:,1),Contours{i}(:,2));
%     %OverlapAll=sum(NumberAll);
%     ContArea=polyarea(Contours{i}(:,1),Contours{i}(:,2));
%     if (ContArea<(2*MaxArea))||(ContArea<(1.5*OverlapAll))
%        PointsInsideActual=sum(inpolygon(EstimatedCellPositions(:,1),EstimatedCellPositions(:,2),Contours{i}(:,1),Contours{i}(:,2))); 
%     if (ContArea<(1.5*OverlapAll))||(PointsInsideActual>0)
%     NumberCluster=inpolygon(ClusterCoord(:,1),ClusterCoord(:,2),Contours{i}(:,1),Contours{i}(:,2));
%     OverlapCluster=sum(NumberCluster);
%     if OverlapAll>1
        if max(size(Contours{i}))>LargeLimit
%             ConvexC=Contours{i}(convhull(Contours{i}(:,1),Contours{i}(:,2)),:);
%             %PIindexList=zeros([];
%             PIindexList=inpolygon(EstimatedCellPositions(:,1),EstimatedCellPositions(:,2),ConvexC(:,1),ConvexC(:,2));
%             
%             PIconvex=sum(PIindexList);
% 
%             if OverlapCluster>5
                LargeContours{p}=Contours{i};
                p=p+1;
%                 PIindexSum=ceil(0.5*(PIindexSum+PIindexList));
%             elseif PIconvex==1
%                 OutsideClusterLarge{r}=Contours{i};
%                 r=r+1;
%             end
        else
            SmallContours{q}=Contours{i};
            q=q+1;
        end
           
end


end

