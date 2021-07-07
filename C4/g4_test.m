
I = imread('test.png');
I = rgb2gray(I);


A = im2double(I);

figure('Name','original image');
imshow(A)


D = dctmtx(size(A,1));
dct = D*A*D';

figure('Name','dct of image');
imshow(dct)




R = D' * dct * D;

R(abs(R)<0.3)=0;

figure('Name','Inverse dct of image');
imshow(R)


