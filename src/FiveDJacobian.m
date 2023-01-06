syms r phi theta1 theta2 theta3
J11 = sin(phi)*sin(theta1)*sin(theta2)*cos(theta3);
J21 = sin(phi)*sin(theta1)*sin(theta2)*sin(theta3);
J31 = sin(phi)*sin(theta1)*cos(theta2);
J41 = sin(phi)*cos(theta1);
J51 = cos(phi);

J12 = r*cos(phi)*sin(theta1)*sin(theta2)*cos(theta3);
J22 = r*cos(phi)*sin(theta1)*sin(theta2)*sin(theta3);
J32 = r*cos(phi)*sin(theta1)*cos(theta2);
J42 = r*cos(phi)*cos(theta1);
J52 = -r*sin(phi);

J13 = r*cos(phi)*cos(theta1)*sin(theta2)*cos(theta3);
J23 = r*cos(phi)*cos(theta1)*sin(theta2)*sin(theta3);
J33 = r*cos(phi)*cos(theta1)*cos(theta2);
J43 = -r*cos(phi)*sin(theta1);
J53 = 0;

J14 = r*cos(phi)*cos(theta1)*cos(theta2)*cos(theta3);
J24 = r*cos(phi)*cos(theta1)*cos(theta2)*sin(theta3);
J34 = -r*cos(phi)*cos(theta1)*sin(theta2);
J44 = 0;
J54 = 0;

J15 = -r*cos(phi)*cos(theta1)*cos(theta2)*sin(theta3);
J25 = r*cos(phi)*cos(theta1)*cos(theta2)*cos(theta3);
J35 = 0;
J45 = 0;
J55 = 0;

J = [J11 J12 J13 J14 J15; J21 J22 J23 J24 J25; J31 J32 J33 J34 J35; J41 J42 J43 J44 J45; J51 J52 J53 J54 J55];
M = det(J)
f = - r^4*cos(phi)^5*cos(theta1)^4*cos(theta2)^3*cos(theta3)^2 - r^4*cos(phi)^5*cos(theta1)^4*cos(theta2)^3*sin(theta3)^2 - r^4*cos(phi)^5*cos(theta1)^4*cos(theta2)*cos(theta3)^2*sin(theta2)^2 - r^4*cos(phi)^5*cos(theta1)^4*cos(theta2)*sin(theta2)^2*sin(theta3)^2 - r^4*cos(phi)^5*cos(theta1)^2*cos(theta2)^3*cos(theta3)^2*sin(theta1)^2 - r^4*cos(phi)^5*cos(theta1)^2*cos(theta2)^3*sin(theta1)^2*sin(theta3)^2 - r^4*cos(phi)^5*cos(theta1)^2*cos(theta2)*cos(theta3)^2*sin(theta1)^2*sin(theta2)^2 - r^4*cos(phi)^5*cos(theta1)^2*cos(theta2)*sin(theta1)^2*sin(theta2)^2*sin(theta3)^2 - r^4*cos(phi)^3*cos(theta1)^4*cos(theta2)^3*cos(theta3)^2*sin(phi)^2 - r^4*cos(phi)^3*cos(theta1)^4*cos(theta2)^3*sin(phi)^2*sin(theta3)^2 - r^4*cos(phi)^3*cos(theta1)^4*cos(theta2)*cos(theta3)^2*sin(phi)^2*sin(theta2)^2 - r^4*cos(phi)^3*cos(theta1)^4*cos(theta2)*sin(phi)^2*sin(theta2)^2*sin(theta3)^2 - r^4*cos(phi)^3*cos(theta1)^2*cos(theta2)^3*cos(theta3)^2*sin(phi)^2*sin(theta1)^2 - r^4*cos(phi)^3*cos(theta1)^2*cos(theta2)^3*sin(phi)^2*sin(theta1)^2*sin(theta3)^2 - r^4*cos(phi)^3*cos(theta1)^2*cos(theta2)*cos(theta3)^2*sin(phi)^2*sin(theta1)^2*sin(theta2)^2 - r^4*cos(phi)^3*cos(theta1)^2*cos(theta2)*sin(phi)^2*sin(theta1)^2*sin(theta2)^2*sin(theta3)^2;
q1 = int(f,r,0,1)
q2 = int(q1,phi,0,pi)
q3 = int(q2,theta1,0,pi)
q4 = int(q3,theta2,0,pi)
q5 = int(q4,theta3,0,2*pi)