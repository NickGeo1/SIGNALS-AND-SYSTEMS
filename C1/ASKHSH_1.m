% Our signal is:
% x(t) = cos(100*pi*t) + cos(200*pi*t) + sin(500*pi*t)
% and our domain for the independent variable (t) is [-10,10].



% 1.Finding the minimum sampling frequency
%
% In order to find the minimum sample frequency (f_{s,min}) , we have to 
% find the maximum frequency implemented in signal.We need to express the 
% given signal as a superposition of fundamental sinusoidal signals.We can 
% write the signal as:
% x(t) = sin(100*pi*t + pi/2) + sin(200*pi*t + pi/2) + sin(500*pi*t)[1]
% 
% In general,a signal may be represented by the
% following mathematical form:
%
% x(t) = A * sin(2*pi*f*t+θ)
%
% A: denotes the amplitude of the oscillation.
% F: denotes the frequency of the oscillation.
% t: denotes the independent time variable.
% Θ: denotes the phase of the oscillation.
% 
% Looking at [1], we can tell that 100*pi = 2*pi*f1, 200*pi = 2*pi*f2 and
% 500*pi = 2*pi*f3.So we have f1=50hz, f2=100hz and f3=250hz.So the max
% frequency implemented on signal is fmax=f3=250hz.
%
% From the Nyquist Sampling Theorem, we know that:
% 
% 1/Ts >= 2*fmax <=> fs >= 2*fmax => fs >= 500hz, where Ts is the sample
% period. So we have f_{s,min} = 500hz





% 2.Representing the signal for -t_max <= t <= t_max where t_max = 10
% and step dt=0.001

dt = 0.001; %step
t_max = 10; %interval bound (-10 <= t <= 10)

t = -t_max:dt:t_max; %t is a 1*(2*t_max/dt + 1)= 1*20.001 dimension matrix that
% holds all the time values. t=[-10,-9.999,...,9.999,10]

x = cos(100*pi*t) + cos(200*pi*t) + sin(500*pi*t); %x is a matrix with
% the same dimensions as t. x holds the signal values, each value
% calculated as the corresponding x(t).
% x=[x(t(1)),x(t(2)),...,x(t(20.001))] = [x(-10),x(-9.999),...,x(10)]
   
figure('Name','exercise 2.2 through 2.5 signal representation'); %makes a figure window
%with the corresponding name

plot(t,x,'-r','LineWidth',1.5); %makes x,y axis representing every value in
%t with its corresponding x(t)

xlabel('-10 \leq t \leq +10'); %x label
ylabel('x(t)'); 
grid on %show grid





% 3.Creating signal with sampling frequency fs=f_{s,min}=500hz <=> Ts=1/f_{s,min} <=>
% Ts=1/500 <=> Ts=0.002s.Original signal can be effeciently reconstructed.
% x(t_max) = x(Nmax*Ts) ==> t_max = Nmax * Ts ==> Nmax = t_max / Ts
% ==> Nmax = 10/0.002 = 5000

Ts = 0.002; %Sampling period
Nmax = t_max / Ts; %Max natural multiple of Ts
n = -Nmax:1:Nmax; %n1=[-Nmax,-Nmax+1,...,Nmax-1,Nmax] (natural multiples of Ts)

xs = cos(100*pi*n*Ts) + cos(200*pi*n*Ts) + sin(500*pi*n*Ts);
% xs is a 1*(2*Nmax + 1) = 1*10001 dimension matrix that contains the 
% signal values for each natural multiple of Ts.
% xs=[xs(n1(1)),xs(n1(2)),...,xs(n1(10001))]=[xs(-5000),xs(-4999),...,xs(5000)]

x1 = zeros(1,length(t)); %x1 is a 1*20.001 dimension matrix with 0.Here,
% we are going to store each reconstructed signal value in the interval of [-tmax,tmax]
% with dt=0.001

% In order to get each reconstructed signal value for each t, we have to calculate the
% sum: 
%        Nmax
%        ----
%        \
% x(t) = /      x(n*Ts)*sinc((t-n*Ts)/Ts), for each t              
%        ----
%       n=-Nmax
%
% Each new x(t) value (let x1(t) be each reconstructed signal value) can be
% calculated with a matrix multiplication, which is equivalent to the
% sum above:
%
% x1(t) = xs*sinc, where 
% xs=[x(-Nmax*Ts),x((-Nmax+1)*Ts),...,x(Nmax*Ts)]
% and 
% sinc = [sinc((t-(-Nmax)*Ts)/Ts),sinc((t-(-Nmax+1)*Ts)/Ts),...,sinc((t-Nmax*Ts)/Ts)]^T
%
%
% Calculating the sum for each x1(t), based on the previusly mentioned 
% matrix multiplication:
for k = 1:1:length(t)
    x1(k) = xs * sinc((t(k)-n*Ts)/Ts)';
end

%representing the results in the same figure window
hold on
plot(t,x1,'-b','LineWidth',1.2);





% 4.Creating signal with sampling frequency fs>f_{s,min}.Let fs=555hz, Ts=1/fs =>
% Ts=1/555 <=> Ts≈0.0018s.Original signal can be effeciently reconstructed.
% x(t_max) = x(Nmax*Ts) ==> t_max = Nmax * Ts ==> Nmax = t_max / Ts
% ==> Nmax = 10/0.0018 ≈ 5556

Ts = 0.0018; %Sampling period
Nmax = t_max / Ts; %Max natural multiple of Ts
n = -Nmax:1:Nmax; %n=[-Nmax,-Nmax+1,...,Nmax-1,Nmax] (natural multiples of Ts)

xs = cos(100*pi*n*Ts) + cos(200*pi*n*Ts) + sin(500*pi*n*Ts);
% xs=[xs(n(1)),xs(n(2)),...,xs(n(11113))]=[xs(-5556),xs(-5555),...,xs(5556)]

x2 = zeros(1,length(t));

%Calculating the sum for each x1(t)
for k = 1:1:length(t)
    x2(k) = xs * sinc((t(k)-n*Ts)/Ts)';
end

%representing the results in the same figure window
hold on
plot(t,x2,'-g','LineWidth',1);





% 5.Creating signal with sampling frequency fs<f_{s,min}.Let fs=125hz, Ts=1/fs =>
% Ts=1/125 <=> Ts=0.008s.Original signal cannot be effeciently reconstructed.
% x(t_max) = x(Nmax*Ts) ==> t_max = Nmax * Ts ==> Nmax = t_max / Ts
% ==> Nmax = 10/0.008 = 1250

Ts = 0.008; %Sampling period
Nmax = t_max / Ts; %Max natural multiple of Ts
n = -Nmax:1:Nmax; %n=[-Nmax,-Nmax+1,...,Nmax-1,Nmax] (natural multiples of Ts)

xs = cos(100*pi*n*Ts) + cos(200*pi*n*Ts) + sin(500*pi*n*Ts);
% xs=[xs(n(1)),xs(n(2)),...,xs(n(2501))]=[xs(-1250),xs(-1249),...,xs(1250)]

x3 = zeros(1,length(t));

%Calculating the sum for each x1(t)
for k = 1:1:length(t)
    x3(k) = xs * sinc((t(k)-n*Ts)/Ts)';
end

%representing the results in the same figure window
hold on
plot(t,x3,'-y','LineWidth',1.2);
ylabel('x(t),x1(t),x2(t) and x3(t)'); %y label




% 6.We observe that reconstructing with the minimum sample frequency
% (fs=500hz Nyquist Sampling Theorem) we can get a quite accurate
% representation of the original signal(blue sample).
% Therefore, reconstructing with a little greater than the minimum
% sample frequency(fs=555hz), we get almost 100% the original signal
% (green sample).Last but not least, we can tell that if we try to
% reconstruct the original signal with sample frequency less than the
% minimum sample frequency(fs=125hz), we get an unwanted result
% (yellow sample).The original signal is being represented with red color.


