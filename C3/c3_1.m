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
T = 0.3; % duration of each note
t=[0:Ts:T];
 
% Signals of all notes
x = sin(2*pi*song*t);

% Add zeros at the end of each note, to create a delay
b=zeros(1,length(song));

delay =200; % number of zeros to be added
 
for i=1:delay
   x(  :,end+1 ) = b;
end

% Reshape x to 1D array
sig = reshape(x',length(song)*(length(t)+delay),1);
% Sends audio signal y to the speaker at sample rate Fs.
soundsc(sig,Fs)
 



