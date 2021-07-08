%1.We consider that the movement between the camera and the captured
%object, lasts T=3s and has constant velocity c=10m/s 

T=3; %movement time
c=10; %velocity
con=3; %a random constant for |K(s)|
s=-2*pi:0.001:2*pi;
K2=con*abs(sin(c*T*s/2)./(c*T*s/2));

%Show |K(s)|
figure('Name','|K(s)|');
plot(s,K2,'-r','LineWidth',1.5);
grid on

%3.


