clc
clear all
close all


% Load image into matrix I.
I = imread('test2.jpg');
% Convert imgage intensity values from uint8 to double.
I = double(I);
% Normalize intensity values from [0..255] to [0..1] interval.
I = I / 255;
% Display original image.
figure('Name','Original Image');
imshow(I);

% Get Red Green and Blue Color Components.
IR = I(:,:,1);
IG = I(:,:,2);
IB = I(:,:,3);




c = 1;
T = 20;

%K = fspecial('motion',c*T,0)
K = 1/(c*T)*ones(1,c*T);

CR = conv2(IR, K);
CG = conv2(IG, K);
CB = conv2(IB, K);


% Compose the original image.
IF(:,:,1) = CR;
IF(:,:,2) = CG;
IF(:,:,3) = CB;

figure('Name','Blurred image');
imshow(IF);
%imwrite(IF, 'mask_image2.jpg');

	% Ensure size is odd so that there is a middle point 
	% about which to rotate the filter by angle.
    
    
    % https://www.quora.com/What-algorithm-is-behind-convolution-in-MATLAB
