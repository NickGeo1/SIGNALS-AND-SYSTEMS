clc
clear all


% Define notes as variables
A = 220;
A_sharp = 220*2^(1/12);
Bb = 220*2^(1/12);
B = 220*2^(2/12);
C = 220*2^(3/12);
C_sharp = 220*2^(4/12);
Db = 220*2^(4/12);
D = 220*2^(5/12);
D_sharp =220*2^(6/12);
Eb = 220*2^(6/12);
E = 220*2^(7/12);
F = 220*2^(8/12);
F_sharp = 220*2^(9/12);
Gb = 220*2^(9/12);
G = 220*2^(10/12);
G_sharp = 220*2^(11/12);
Ab = 220*2^(11/12);
 

% This is a sample song
song = [A;A;E;E;F_sharp;F_sharp;E;E;D;D;C_sharp;C_sharp;B;B;A;A];



Fs=8000; % Sampling frequency
Ts=1/Fs;

sig=[]
A=1 % volume of note
for i = 1:numel(song)
    T=0.3;
    t=[0:Ts:T];
    A = A*0.7; % volume reduced to 70%
    x = A*sin(2*pi*song(i)*2*t);


    delay =200; % number of zeros to be added

    for y=1:delay
       x(  :,end+1 ) = 0;
    end
    

    for j=1:numel(x)
        sig = [sig;x(j)];%append elements of x to sig
    end


end

%play the signal
sig(1)
soundsc(sig,Fs)





