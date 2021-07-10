% Clear all memory and command line
% Close all figures
clc
clear all
close all

% Load a builtin demo image from matlab into the matrix I
I = imread('peppers.png');
% Convert imgage intensity values from uint8 to double.
I = double(I);
% Normalize values of I from [0..255] to [0..1] interval.
I = I / 255;
% Get Red Green and Blue Color Components.
IR = I(:,:,1);
IG = I(:,:,2);
IB = I(:,:,3);

% The Discrete Cosine Transform (DCT) of a 2-d signal g(i, k) with 0<=i
% and k<=N-1, can be expressed with orthonormal matrixes:
%
%           G_c = C' * g * C        (1)
% 
% where C(i, m) = a(m)*cos((pi*(2*i + 1)*m)/(2*N))
%
%               = a(m)*cos(((i+0.5)*pi*m)/N)
%
% and a(0)=1/sqrt(N), a(m) = sqrt(2/N).
%
% The inverse DCT of a 2-d signal is given by:
%
%           g = C * G_c * C'        (2)
%

% To compress an image, first compute the DCT (1) of an 8x8 block from the
% image. Discard some of the DCT coefficients and then reconstruct
% the image using the  inverse DCT (2). This is done for every 8x8 block of
% the image. The coefficients choosen in each block are defined by another 8x8 array called 'mask'.
% The more 'zeros' the mask has, the more compressed the image will be.

% This is the principle of the JPEG image compression algorithm. The input
% image acn also be divided by blocks of 16x16.

% Calculate C of the DCT formula
N = 8;
C = zeros(N);

for i = 0:(N-1)
    for m = 0:(N-1)
        if m == 0
            a = sqrt(1/N);
        else
            a = sqrt(2/N);
        end
        C(i+1,m+1) = a*cos((i+0.5)*pi*m/N);
    end
end

% Define the function that will be used for each block.
% In this case the function is formula (1)
dct = @(block_struct) C' * block_struct.data * C;

% Use blockproc function to calculate the result of dct for each 8x8 block
% of IR/IG/IB. The results pass to DCTR/DCTG/DCTB accordingly.
DCTR = blockproc(IR,[8 8],dct);
DCTG = blockproc(IG,[8 8],dct);
DCTB = blockproc(IB,[8 8],dct);


% Define four masks with a different number of ones.

mask1 = [1   1   1   1   1   0   0   0
         1   1   1   1   0   0   0   0
         1   1   1   0   0   0   0   0
         1   1   0   0   0   0   0   0
         1   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0];
    
    


mask2 = [1   1   1   1   0   0   0   0
         1   1   1   0   0   0   0   0
         1   1   0   0   0   0   0   0
         1   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0];
    
    

mask3 = [1   1   1   0   0   0   0   0
         1   1   0   0   0   0   0   0
         1   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0];
    
     

mask4 = [1   1   0   0   0   0   0   0
         1   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0
         0   0   0   0   0   0   0   0];
     

% Define the function that will be used for each block.
% In this case the function is formula (2)
invdct = @(block_struct) C * block_struct.data * C'; 


%%%%%%%%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%% 
% First, for each 8x8 block of DCTR/DCTG/DCTB mutliply by value with mask1 
% to keep the specific DCT coefficients 
DCTR1 = blockproc(DCTR,[8 8],@(block_struct) mask1 .* block_struct.data);
DCTG1 = blockproc(DCTG,[8 8],@(block_struct) mask1 .* block_struct.data);
DCTB1 = blockproc(DCTB,[8 8],@(block_struct) mask1 .* block_struct.data);

% Then calculate the invdct of each block by the defined function
IDCTR1 = blockproc(DCTR1,[8 8],invdct);
IDCTG1 = blockproc(DCTG1,[8 8],invdct);
IDCTB1 = blockproc(DCTB1,[8 8],invdct);

% Compose the output image.
IDCT1(:,:,1) = IDCTR1;
IDCT1(:,:,2) = IDCTG1;
IDCT1(:,:,3) = IDCTB1;
%%%%%%%%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%% 

% REPEAT PROSSES FOR ALL THE MASKS

%%%%%%%%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%          
DCTR2 = blockproc(DCTR,[8 8],@(block_struct) mask2 .* block_struct.data);
DCTG2 = blockproc(DCTG,[8 8],@(block_struct) mask2 .* block_struct.data);
DCTB2 = blockproc(DCTB,[8 8],@(block_struct) mask2 .* block_struct.data);

IDCTR2 = blockproc(DCTR2,[8 8],invdct);
IDCTG2 = blockproc(DCTG2,[8 8],invdct);
IDCTB2 = blockproc(DCTB2,[8 8],invdct);

% Compose the output image.
IDCT2(:,:,1) = IDCTR2;
IDCT2(:,:,2) = IDCTG2;
IDCT2(:,:,3) = IDCTB2;
%%%%%%%%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%% 


 
%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  
DCTR3 = blockproc(DCTR,[8 8],@(block_struct) mask3 .* block_struct.data);
DCTG3 = blockproc(DCTG,[8 8],@(block_struct) mask3 .* block_struct.data);
DCTB3 = blockproc(DCTB,[8 8],@(block_struct) mask3 .* block_struct.data);

IDCTR3 = blockproc(DCTR3,[8 8],invdct);
IDCTG3 = blockproc(DCTG3,[8 8],invdct);
IDCTB3 = blockproc(DCTB3,[8 8],invdct);

% Compose the output image.
IDCT3(:,:,1) = IDCTR3;
IDCT3(:,:,2) = IDCTG3;
IDCT3(:,:,3) = IDCTB3;
%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%



%%%%%%%%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  
DCTR4 = blockproc(DCTR,[8 8],@(block_struct) mask4 .* block_struct.data);
DCTG4 = blockproc(DCTG,[8 8],@(block_struct) mask4 .* block_struct.data);
DCTB4 = blockproc(DCTB,[8 8],@(block_struct) mask4 .* block_struct.data);

IDCTR4 = blockproc(DCTR4,[8 8],invdct);
IDCTG4 = blockproc(DCTG4,[8 8],invdct);
IDCTB4 = blockproc(DCTB4,[8 8],invdct);

% Compose the output image.
IDCT4(:,:,1) = IDCTR4;
IDCT4(:,:,2) = IDCTG4;
IDCT4(:,:,3) = IDCTB4;
%%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%%  %%%%%%% 



% Display all matrixes
figure('Name','Original Image');
imshow(I)

figure('Name','DCT compressed image with mask1');
imshow(IDCT1)

figure('Name','DCT compressed image with mask2');
imshow(IDCT2)

figure('Name','DCT compressed image with mask3');
imshow(IDCT3)

figure('Name','DCT compressed image with mask4');
imshow(IDCT4)

% Save images to check file size.               SIZE

imwrite(I, 'peppers.png');                      %281K
imwrite(IDCT1, 'compressed_peppers_mask1.png'); %243K
imwrite(IDCT2, 'compressed_peppers_mask2.png'); %228K
imwrite(IDCT3, 'compressed_peppers_mask3.png'); %208K
imwrite(IDCT4, 'compressed_peppers_mask4.png'); %172k
