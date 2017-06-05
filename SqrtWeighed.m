function [WeighedImage] = SqrtWeighed(Image)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if max(max(Image))<=1
    Image=round(255*double(Image));
end
Image=double(Image);
WeighedImage=round(255.*(Image.^0.5)./(255.^0.5));
end
