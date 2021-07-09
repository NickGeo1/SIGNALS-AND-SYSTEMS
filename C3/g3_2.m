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

sig1=[] % signal to be heard, will now be calculated with a loop, is one octave up


% for each note of the song, duration is changed
for i = 1:numel(song)
    T=0.1*i; % current duration
    t=[0:Ts:T];
    x = sin(4*pi*song(i)*t); % mutiply by 2 for one octave up


    delay =200; % number of zeros to be added

    for y=1:delay
       x(  :,end+1 ) = 0;
    end
    

    for j=1:numel(x)
        sig1 = [sig1;x(j)];%append elements of x to sig
    end


end

%play the signal
sig1(1)
soundsc(sig1,Fs)


pause(numel(sig1)/8000+1) %sig1 plays for number_of_samples / 8000 seconds
 

sig2=[]% signal to be heard, is one octave down
% similar loop to the first one
for i = 1:numel(song)
    T=0.1*i;
    t=[0:Ts:T];
    x = sin(pi*song(i)*t); %divide by 2 to go one octave down


    delay =200; % number of zeros to be added

    for y=1:delay
       x(  :,end+1 ) = 0;
    end
    

    for j=1:numel(x)
        sig2 = [sig2;x(j)];
    end


end

%play the signal
sig2(1)
soundsc(sig2,Fs)


