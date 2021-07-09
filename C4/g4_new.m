clear all
I = imread('cameraman.tif');
I = double(I);
I = I / 255;



A = zeros(8);


for i = 0:7
    for j = 0:7
        if i == 0
            a = sqrt(1/8);
        else
            a = sqrt(2/8);
        end
        A(i+1,j+1) = a*cos((j+0.5)*pi*i/8);
    end
end


T = A;
dct = @(block_struct) T * block_struct.data * T';
B = blockproc(I,[8 8],dct);




mask = [1   1   1   0   0   0   0   0
        1   1   0   0   0   0   0   0
        1   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0
        0   0   0   0   0   0   0   0];
    
    

        
        
B2 = blockproc(B,[8 8],@(block_struct) mask .* block_struct.data);



invdct = @(block_struct) T' * block_struct.data * T;
I2 = blockproc(B2,[8 8],invdct);

imshow(I)
figure
imshow(I2)




