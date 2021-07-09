clc
clear all
close all

%1.We consider that the movement between the camera and the captured
%object, lasts T=3s and has constant velocity c=10m/s 

T=3; %movement time
c=10; %velocity
con=3; %a random constant for |K(s)|
s=-pi:0.001:pi; %signal construction sample frequency = 0.001.Interval [-pi,pi] 
K=con*abs(sin(c*T*s/2)./(c*T*s/2)); %calculate |K(s)| for every s

%Show |K(s)| with c=10
figure('Name','|K(s)|');
plot(s,K,'-r','LineWidth',1.5);
grid on

%Show |K(s)| with c=20
c=20;
K=con*abs(sin(c*T*s/2)./(c*T*s/2));
hold on
plot(s,K,'-b','LineWidth',1);
xlabel('Frequency s');
ylabel('F(s) with c=10 (red) and with c=20 (blue)')





% 2.Initialy, we will solve for s the equation |K(s)|=0 (s!=0) and we will
% keep the first positive solution:

% |K(s)|=0 <=> sin(c*T*s/2) = 0
% <=> 
% c*T*s/2 = 2*k*pi      |     | s = 4*k*pi/c*T
% or                    | <=> | or                       k ∈ Z* 
% c*T*s/2 = 2*k*pi + pi |     | s = 2*pi*(2*k + 1)/c*T
%  
% for k = 0 we have: s0=0 (solution denied) or s0 = 2*pi/c*T (accepted)
% 
% We observe that as greater as the c is, so smaller the s0 becomes, and
% vice versa. That means, that as the c becames greater, we have a lack
% of low frequencies in the output image(smaller frequencies make K(s)=0).
% So, as we increase the velocity, the blur also increases.





%3.
% Load image into matrix I.
img = imread('tree.png');
% Convert imgage intensity values from uint8 to double.
img = double(img);
% Normalize intensity values from [0..255] to [0..1] interval.
img = img / 255;
% Display original image.
figure('Name','Original Image');
imshow(img);

% Get Red Green and Blue Color Components.
imgR = img(:,:,1);
imgG = img(:,:,2);
imgB = img(:,:,3);

c = 10; %velocity
T = 3; %movement time

K = 1/(c*T)*ones(1,floor(c*T)); %K is a continuous function, therefore, we have 
% to take some discrete values of it in order to perform the following
% convolution. We take ⌊c*T⌋ values, each value equal to 1/(c*T), because function
% is constant and equal to 1/(c*T) for 0<=x<=c*T.

CR = conv2(imgR, K);
CG = conv2(imgG, K);
CB = conv2(imgB, K);%Convolution between one two-dimensional discrete time signal(image),
% and one-dimensional discrete time signal(K values).Each convolution is
% calculated based on the following formula:
%
%            p    p
%           ---  ---
%           \    \
% C(m,n) =  /    /    = A(i,j)*B(m−i+1,n−j+1)
%           ---  ---
%           i=1  j=1
%
%
% source:
% https://www.quora.com/What-algorithm-is-behind-convolution-in-MATLAB
%
%
% Compose the original image.
blured_img(:,:,1) = CR;
blured_img(:,:,2) = CG;
blured_img(:,:,3) = CB;

figure('Name','Blurred image');
imshow(blured_img);%show blured img   
    

