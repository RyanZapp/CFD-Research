clear
close all
clc
% Variables that can be changed
Chord_Lengths = [0.0762, 0.0762, 3]; % chord lengths in inches
bladeName = 'NREL_5MW'; % Change this to match your blade name
numberBlades = '3';
H = 2.5; % Half the height of the shaft (non-varying parameter)
N1 = 34; % Number of airfoil segments used to approximate the blade
S = 0.5; % Value that you scale TAxis by
z = linspace(0,2*H,N1)';
AR = 1.5; % Aspect Ratio...This value is eqaul to H/R;
R_true = H/AR;
% End of Variables that can be changed

space = ' ';
z1 = linspace(-H,H,N1)';

N2 = 1000;
aa = linspace(1,5,N2);
R_approx = -aa + aa.*cosh(H./aa);
R = interp1(R_approx,R_approx,R_true,'nearest');
for j = 1:N2
    if (R_approx(j) == R)
        a = aa(j);
    end
end
% The block below uses a search to find the value of a
% for j = 1:N2
%     if (R_approx(j) <= 1.67)&(R_approx(j) >= 1.6659)
%         disp('Potential Value of a:');
%         disp(aa(j));
%         a = aa(j);
%     end
% end
% The equation below is for visualization purposes
% The equation we will export must be shifted and flipped due to the way
% that Qblade takes in data
y_viz = a*cosh(z1/a) - a*cosh(H/a);
y_calc = -a*cosh(z/a - H/a) + a*cosh(H/a); % use this for calculations in Qblade
figure
plot(z1,y_viz,'b','LineWidth',1.5)
hold on
plot([-H H],[0 0],'k','LineWidth',1.5)
hold off
title('Catenary Blade Shape Plot (For Visualization Purposes)')
grid on
axis equal
axis([-(H+1) (H+1) -inf 1])
%maximum value occurs -a + a cosh(H/a)
Twist = num2str(zeros(N1,1));
Circ = num2str(zeros(N1,1));
Height = num2str(0.3048*z);
Radius = num2str(0.3048*y_calc);
TAxis = num2str(S*ones(N1,1));

Chord = ones(N1,1);

Chord(1:5) = Chord_Lengths(1)*Chord(1:5);
Chord(6:29) = Chord_Lengths(2)*Chord(6:29);
Chord(30:34) = Chord_Lengths(1)*Chord(30:34);
%Chord(5:10) = Chord_Lengths(2)*Chord(5:10);
%Chord(11:15) = Chord_Lengths(3)*Chord(11:15);

% In between the line above and the line below, make assignments like
% the one made in the line above in order to create the chord lengths you
% want
% I reallty only need to go halfway through the chord..i.e. to N1/2
% Then I can just reverse everything on the way back
Chord = num2str(Chord);
%for i = 1:N1
%T(i,:) = {Height(i),Chord(i),Twist(i),1,1,TAxis(i),'dog'};
%end
T1 = 'HEIGHT [m]'; % This used to read "POS [m]"
T2 = 'CHORD [m]';
T3 = 'RADIUS [m]'; % delete these if this idea does not work
T4 = 'TWIST [deg]';
T5 = 'CIRC [deg]';
T6 = 'TAXIS [-]';
%T7 = 'POLAR_FILE';

s1 = '----------------------------------------QBlade Blade Definition File------------------------------------------------';
s2 = '----------------------------------------Object Name-----------------------------------------------------------------';
s3 = '----------------------------------------Parameters------------------------------------------------------------------';
s4 = '----------------------------------------Blade Data------------------------------------------------------------------';

k = zeros(2,1);
timeSpec = 'Time : %s\n';
Time = datetime('now','TimeZone','local','Format','HH:mm:ss');
dateSpec = 'Date : %s\n';
Date = datetime('now','TimeZone','local','Format','dd.MM.yyyy');
%formatSpec = '%.5f';

% Name your Blade File Below
fid = fopen('BladeFile2.txt','w');

fprintf(fid,s1);
fprintf(fid,'\n');
fprintf(fid,'Generated with : QBlade CE v 2.0.4.8_alpha windows\n');
fprintf(fid,'Archive Format: 310002\n');
fprintf(fid,timeSpec,Time);
fprintf(fid,dateSpec,Date);
fprintf(fid,'\n');
fprintf(fid,s2);
fprintf(fid,'\n');
fprintf(fid,bladeName);
for j = 1:(41 - length(bladeName))
    fprintf(fid,space);
end
fprintf(fid,'OBJECTNAME');
for j = 1:9
    fprintf(fid,space);
end
fprintf('- the name of the blade object\n');
fprintf(fid,'\n');
fprintf(fid,s3);
fprintf(fid,'\n');
fprintf(fid,'VAWT');
for j = 1:37
    fprintf(fid,space);
end
fprintf(fid,'ROTORTYPE');
for j = 1:10
    fprintf(fid,space);
end
fprintf(fid,'- the rotor type\n');
fprintf(fid,numberBlades);
for j = 1:40
    fprintf(fid,space);
end
fprintf(fid,'NUMBLADES');
for j = 1:10
    fprintf(fid,space);
end
fprintf(fid,'- number of blades\n');
fprintf(fid,s4);
fprintf(fid,'\n');
fprintf(fid,T1);
for j = 1:(20-length(T1))
    fprintf(fid,space);
end
fprintf(fid,T2);
% spacings of 20
for j = 1:(20-length(T2))
    fprintf(fid,space);
end
fprintf(fid,T3);
for j = 1:(20-length(T3))
    fprintf(fid,space);
end
fprintf(fid,T4);
for j = 1:(20-length(T4))
    fprintf(fid,space);
end
fprintf(fid,T5);
for j = 1:(20-length(T5))
    fprintf(fid,space);
end
fprintf(fid,T6);
for j = 1:(20-length(T6))
    fprintf(fid,space);
end
%fprintf(fid,T7);
fprintf(fid,'\n');
for j = 1:N1
    fprintf(fid,Height(j,:));
    for i = 1:(20-length(Height(j,:)))
        fprintf(fid,space);
    end
    fprintf(fid,Chord(j,:));
    for i = 1:(20-length(Chord(j,:)))
        fprintf(fid,space);
    end
    fprintf(fid,Radius(j,:));
    for i = 1:(20-length(Radius(j,:)))
        fprintf(fid,space);
    end
    fprintf(fid,Twist(j,:));
    for i = 1:(20-length(Twist(j,:)))
        fprintf(fid,space);
    end
    fprintf(fid,Circ(j,:));
    for i = 1:(20-length(Circ(j,:)))
        fprintf(fid,space);
    end
    fprintf(fid,TAxis(j,:));
    for i = 1:(20-length(TAxis(j,:)))
        fprintf(fid,space);
    end
    %fprintf(fid,'s'); % Replace with polar file later (also, turn it into a for loop)
    fprintf(fid,'\n');
end
% I have removed the need for a polar file, but i can add it back in if it
% becomes necessary