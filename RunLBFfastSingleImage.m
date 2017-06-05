function [CellChs,HandledImages] = RunLBFfastSingleImage(Channels,Clusters,SigmaL,iterNum,FigureHandle)
                                             % (Channels,Clusters,SigmaL,iterNum,FigureHandle,Mode,Background);
% This Matlab file demomstrates a level set method in Chunming Li et al's paper
%    "Minimization of Region-Scalable Fitting Energy for Image Segmentation",
%    IEEE Trans. Image Processing, vol. 17 (10), pp.1940-1949, 2008.
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://www.engr.uconn.edu/~cmli/
%
% Note 1: The original model (LBF) with a small scale parameter sigma, such as sigma = 3, is sensitive to 
%         the initialization of the level set function. Appropriate initial level set functions are given in 
%         this code for different test images. 
% Note 2: There are several ways to improve the original LBF model to make it robust to initialization.
%         One of the improved LBF algorithms is implemented by the code in the following link:
%                http://www.engr.uconn.edu/~cmli/code/LBF_v0.1.rar


    Img=Channels;
    ClusterMask=Clusters-0.5;
    initialLSF=ClusterMask;
    HandledImages=uint8(Img);
    %toc
%      disp('pos 2')
    CellChs=RunLevelSetFast(Img,initialLSF,SigmaL,iterNum,FigureHandle,(i-1));
%     disp('pos 3')


end

function [Cells]=RunLevelSetFast(Img,initialLSF,Sigma0,iterNum,FigureHandle,ChNo)
    
%tic
% disp('pos 3')
    lambda1 = 1.0;  lambda2 = 1.0; nu = 0.001*255*255;
   
    c0 = 2;
    u=c0*sign(initialLSF);    
    Img=double(Img(:,:,1));
    %axes(FigureHandle); % Switches focus to this axes object.
    popupfigure = figure;
    imagesc(Img, [0, 255]);colormap(gray);hold on; axis off;
    [nrow, ncol]=size(Img);
    contour(u,[0 0],'r');
    title(['Channel ' mat2str(ChNo) ' Cell Search Initial State'])
    pause(.05);
    pause(0.1);
    timestep = .1;
    mu = 1;
    epsilon = 1.0;
    sigma=Sigma0;
    %K=fspecial('gaussian',round(2*sigma)*2+1,sigma); % Gaussian kernel
    %KI=conv2(Img,K,'same');
    %KONE=conv2(ones(nrow, ncol),K,'same');

%timestep = .1;% time step
%mu = 1;% coefficient of the level set (distance) regularization term P(\phi)

%epsilon = 1.0;% the papramater in the definition of smoothed Dirac function
%sigma=3.0;    % scale parameter in Gaussian kernel
              % Note: A larger scale parameter sigma, such as sigma=10, would make the LBF algorithm more robust 
              %       to initialization, but the segmentation result may not be as accurate as using
              %       a small sigma when there is severe intensity inhomogeneity in the image. If the intensity
              %       inhomogeneity is not severe, a relatively larger sigma can be used to increase the robustness of the LBF
              %       algorithm.
 %sz=[round(2*sigma)*2+1;round(2*sigma)*2+1];
              %sz = sz(:)';  % force column vector
 
 %szHalf = (sz-1)/2;
%if (~isequal(szHalf,round(szHalf)))
 %   sz=[(sigma+1);(sigma+1)];
 %   szHalf = (sz-1)/2;
%end
              
K=fspecial('gaussian',round(2*sigma)*2+1,sigma);     % the Gaussian kernel
%K{2}=gausskernel(szHalf(2),sigma);
%Ksigma=
I = Img;
KI=convolve2(Img,K,'same');     % compute the convolution of the image with the Gaussian kernel outside the iteration
                            % See Section IV-A in the above IEEE TIP paper for implementation.
%KI=convnsepB(K{:}, padreplicate(Img,(sz-1)/2), 'valid');                                              
KONE=convolve2(ones(size(Img)),K,'same');  % compute the convolution of Gassian kernel and constant 1 outside the iteration
                                       % See Section IV-A in the above IEEE TIP paper for implementation.
%KONE=convnsepB(K{:}, padreplicate(ones(size(Img)),(sz-1)/2), 'valid');  
% start level set evolution
%toc
 %disp('pos 4')
for n=1:iterNum
    u=LSE_LBF_fastI(u,I,K,KI,KONE, nu,timestep,mu,lambda1,lambda2,epsilon,1);
    if mod(n,20)==0
        pause(0.1);
        imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
        [c,h] = contour(u,[0 0],'r');
        iterNum=[num2str(n), ' iterations'];
        title(iterNum);
        hold off;
    end
    %toc %this is the toc that takes time
end
imagesc(Img, [0, 255]);colormap(gray);hold on;axis off,axis equal
[c,h] = contour(u,[0 0],'r');
totalIterNum=[num2str(n), ' iterations'];
title(['Final contour, ', totalIterNum]);

Cells=u;
%toc
% disp('pos 5')
%figure;
%mesh(u);
%title('Final level set function');

close(popupfigure);

end

function u = LSE_LBF_fastI(u0,Img,Ksigma,KI,KONE,nu,timestep,mu,lambda1,lambda2,epsilon,numIter)
% LSE_LBF implements the level set evolution (LSE) for the method in Chunming Li et al's paper:
%       "Minimization of Region-Scalable Fitting Energy for Image Segmentation", 
%        IEEE Trans. Image Processing(TIP), vol. 17 (10), pp.1940-1949, 2008.
%
% Author: Chunming Li, all rights reserved
% E-mail: li_chunming@hotmail.com
% URL:  http://www.engr.uconn.edu/~cmli/
%
% For easy understanding of my code, please read the comments in the code that refer
% to the corresponding equations in the above IEEE TIP paper. 
% (Comments added by Ren Zhao at Univ. of Waterloo)

u=u0;
for k1=1:numIter
    u=NeumannBoundCond(u);
    K=curvature_central(u);                             

    DrcU=(epsilon/pi)./(epsilon^2.+u.^2);               % eq.(9)

    [f1, f2] = localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon);


    %%% compute lambda1*e1-lambda2*e2
    s1=lambda1.*f1.^2-lambda2.*f2.^2;                   % compute lambda1*e1-lambda2*e2 in the 1st term in eq. (15) in IEEE TIP 08
    s2=lambda1.*f1-lambda2.*f2;
    dataForce=(lambda1-lambda2)*KONE.*Img.*Img+convolve2(s1,Ksigma,'same')-2.*Img.*convolve2(s2,Ksigma,'same');
   %  Sigma = Sigma(:)';  % force column vector
 
 %   szHalf = (Sigma-1)/2;
    
    %     kernel{1} = gausskernel(szHalf(1),arg);
    %     kernel{2} = gausskernel(szHalf(2),arg);
    %convnsepB(kernel{:}, padreplicate(data,(sz-1)/2), 'valid');
                                                        % eq.(15)
    A=-DrcU.*dataForce;                                 % 1st term in eq. (15)
    P=mu*(4*del2(u)-K);                                 % 3rd term in eq. (15), where 4*del2(u) computes the laplacian (d^2u/dx^2 + d^2u/dy^2)
    L=nu.*DrcU.*K;                                      % 2nd term in eq. (15)
    u=u+timestep*(L+P+A);                               % eq.(15)
end
end

function [f1, f2]= localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon)
% compute f1 and f2
Hu=0.5*(1+(2/pi)*atan(u./epsilon));                     % eq.(8)

I=Img.*Hu;
c1=convolve2(Hu,Ksigma,'same'); % Hu,'same');                             
c2=convolve2(I,Ksigma,'same'); %I,'same');                              % the numerator of eq.(14) for i = 1
f1=c2./(c1);                                            % compute f1 according to eq.(14) for i = 1
f2=(KI-c2)./(KONE-c1);                                  % compute f2 according to the formula in Section IV-A, 
end                                                        % which is an equivalent expression of eq.(14) for i = 2.
                                                         

function g = NeumannBoundCond(f)
% Neumann boundary condition
[nrow,ncol] = size(f);
g = f;
g([1 nrow],[1 ncol]) = g([3 nrow-2],[3 ncol-2]);  
g([1 nrow],2:end-1) = g([3 nrow-2],2:end-1);          
g(2:end-1,[1 ncol]) = g(2:end-1,[3 ncol-2]);  
end

function k = curvature_central(u)                       
% compute curvature
[ux,uy] = gradient(u);                                  
normDu = sqrt(ux.^2+uy.^2+1e-10);                       % the norm of the gradient plus a small possitive number 
                                                        % to avoid division by zero in the following computation.
Nx = ux./normDu;                                       
Ny = uy./normDu;
[nxx,~] = gradient(Nx);                              
[~,nyy] = gradient(Ny);                              
k = nxx+nyy;                                            % compute divergence
end

function b=padreplicate(a, padSize)
numDims = length(padSize);
idx = cell(numDims,1);
 for k = 1:numDims
 M = size(a,k);
   onesVector = ones(1,padSize(k));
   idx{k} = [onesVector 1:M M*onesVector];
 end
 b = a(idx{:});
end