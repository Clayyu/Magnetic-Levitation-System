% This program is designed to simulate data used to calculate the "alpha" 
% and "beta" in order to determine the transfer function of the MLS System
% as well as the percent overshoot and settling time of the PID Controller.

% d = Distance in millimeters
d_mm = [19.00
19.50
20.00
20.50
21.00
21.50
22.00
22.50
23.00
23.50
24.00
24.50
25.00
25.50
26.00
26.50
27.00
27.50
28.00
28.50
29.00
29.50
30.00
30.50
31.00
31.50
32.00
32.50
33.00
33.50
34.00
34.50
35.00
35.50
36.00
36.50
37.00
37.50
38.00
38.50
39.00
39.50
40.00
40.50
41.00
41.50
42.00
42.50
43.00
43.50
44.00
44.50
45.00
45.50
46.00
46.50
47.00
47.50
48.00
48.50
49.00
49.50
50.00
50.50
51.00
51.50
52.00
52.50
53.00
53.50
54.00
54.50
55.00
55.50
56.00
56.50
57.00
57.50
58.00
58.50
59.00
59.50
60.00
60.50
61.00
61.50
62.00
62.50
63.00
63.50
64.00
64.50
65.00
65.50
66.00
66.50
67.00
67.50
68.00
68.50
69.00
69.50
70.00
70.50
71.00
71.50
72.00
72.50
73.00
73.50
];

% V = Voltage measured in volts
V = [0.0740
0.1850
0.1900
0.1630
0.2000
0.3280
0.3340
0.3000
0.3700
0.3800
0.4070
0.4100
0.4350
0.4580
0.4760
0.5100
0.5450
0.5650
0.5960
0.6260
0.6480
0.6790
0.7140
0.7460
0.7580
0.7840
0.8130
0.8540
0.8700
0.9000
0.9270
0.9540
0.9890
1.0220
1.0520
1.0880
1.1220
1.1510
1.2090
1.2260
1.2610
1.3100
1.3400
1.3900
1.4270
1.4180
1.4470
1.4720
1.5250
1.5800
1.6540
1.7400
1.7490
1.8000
1.8900
1.9660
2.0400
2.1150
2.2200
2.3280
2.5000
2.6400
2.5300
2.6950
2.8450
2.9500
3.0000
3.3100
3.5740
3.7230
3.7800
4.0500
4.7300
4.5950
4.8450
5.0730
5.4200
5.5650
5.8150
6.1100
6.4600
6.6900
6.8700
7.3000
7.7100
8.0900
8.4300
8.7300
9.1900
9.5800
9.9300
10.4000
10.9200
11.3750
11.6800
12.3000
12.8250
13.3900
13.7600
14.4100
15.1400
15.6200
16.1900
16.9000
17.6900
18.4200
18.8900
19.8000
20.6600
21.5200
];

% Resistance of coil in Ohms before measurements were taken
R = 3.90;

% Current of the coil calculated in Amperes
Ic = V/R;

% d_M = Distance converted from millimeters to meters
d_M = d_mm*.001;

% m2a = Mass of ball in kilograms
m2a = 0.016;

% g = Gravity in m/s^2
g = 9.81;

% y = polynomial equation of relationship between current in the coil and
% weight of the ball
y = (Ic.^2)/(m2a*g);

% p5 = Polyfit used to morph data based on nth polynomial degree
p5 = polyfit(d_M,y,5);

% Ii = Symbolic variable used to differentiate with respect to the Ii
% z = Symbolic variable used to express p5s in terms of the z
% z = Symbolic variable used to differentiate with respect to the z
% z = sym('z')
syms Ii z;

% p5s = Polyfit variable converted to symbolic data type
p5s = poly2sym(p5,z);

% p5ss = Polyfit variable evaluated numerically using the vpa function to
% approximate fractions.
p5ss = vpa(p5s);

% fiz = function f(i,z) of equation y.
fiz = g-(Ii^2/(m2a*p5ss));

% differentiate function f(i,z) with respect to the z.
dfdz = diff(fiz,z);

% differentiate function f(i,z) with respect to the Ii.
dfdi = diff(fiz,Ii);

% dfdzv = dfdz evaluated numerically using the vpa function to approximate
% fractions
dfdzv = vpa(dfdz);

% dfdiv = dfdi evaluated numerically using the vpa function to approximate
% fractions
dfdiv = vpa(dfdi);

% Ii = 2.162 = The current in the coil at 6.2 centimeters expressed in
% amperes.
Ii = 2.162;

% z = 0.062 = The distance from the coil at 2.162 amperes expressed in
% meters.
z = 0.062;

% alpha = 1.8951e+03 = alpha of transfer function for the MLS system.
alpha = eval(dfdzv)

% beta = -9.4666 = beta of transfer function for the MLS system.
beta = eval(dfdiv)

% MLStf = -9.467 / (s^2-1895) = transfer function of the MLS system.
MLStf = tf(beta,[1 0 -alpha])

% MLStfpzmap = 43.5329 and -43.5329 = the poles of the MLS system.
MLStfpzmap = pzmap(MLStf)

% Root Locus of the "Uncompensated MLS System"
rlocus(MLStf),title('Root Locus of the "Uncompensated MLS System"')

pause

% wn^2/(s^2+(2*zeta*wn*s)+wn^2) = The equation taken into consideration for
% designing the stability of the system.

% The following values for Percent Overshoot and the Settling Time are 
% chosen to determine the stability of the system.

% PO = Percent Overshoot = 3mm/62mm = 5%
PO = 5;

% Note: 3 mm overshoot was chosen because 2 mm overshoot was considered
% too small of an overshoot that could be achieved with the stability
% naturally inherent in this type of system.
% It is believed that this 3mm overshoot is within acceptable tolerance for
% the system.

% ts = Settling Time
ts = 1.5

zx = (log(PO/100))^2;

% zeta = 0.6901
zeta = sqrt(zx/(pi^2+zx))

% wn = 3.8641
wn = 4/(zeta*ts)

% s = -2.6667+2.7965i = pole produced from system criteria(ts, PO, etc..)
s = -zeta*wn+j*wn*sqrt(1-zeta^2)

% g = 9.467 / (s^2-1895) = transfer function of the MLS system.
g = tf(-beta,[1 0 -alpha])

% h = feedback of the MLS system without the hall effect sensor is unity
% gain or one.
h = tf([1],[1]);

% pic = (s+4)/s
% pic = transfer function of the proportional integral used to develop the
% PID controller.
pic = tf([1 4],[1 0])

% sysu = (9.467s+37.87) / (s^3-1895s)
% sysu = transfer function of the MLS system with the PIC.
sysu = g*h*pic;

% Root Locus of the "Uncompensated System with the PI"
rlocus(sysu),title('Root Locus of the Uncompensated System with the PI')

pause

% p = -0.0014+0.0038i = angle of pole based on system criteria
p = evalfr(sysu,s)

% ap = 110.4196 = pole angle converted to degrees
ap = (180/pi)*angle(p)

% apt = 69.5804 = Total angle needed to put pole on root locus
apt = 180-ap

% x = 1.0411 = Distance on real axis that puts pole on root locus
% that then produces the zero location for the proportional-derivative.
x=(imag(s)/tan((pi/180)*apt))

% pd = s+3.708 = transfer function of the proportional derivative used to
% develop the PID controller.
pd = tf([1 -real(s)+x],1)

% rlocus(sysu*pd) = Root Locus of the "Compensated System with the
% Proportional Integral and the Proportional Derivative."
rlocus(sysu*pd),title('Root Locus of the "Compensated System with the Proportional Integral and the Proportional Derivative.')

pause

% [kd, pp] = function used to generate a kd, or gain, value for the PID
% controller from the root locus compensated with the PD and PI.
[kd, pp] = rlocfind(sysu*pd,s);

% kd = 83.7088
kd

% mult = 2 = The Variable Defined to Tune the Kd, or Gain, for the PID
% Controller.
mult = 2;

% kdmult = 167.4176 = The Total Kd, or gain, of the PID Controller.
kdmult = kd*mult

% sysPID = (1585s^2+1.222e04s+2.351e04) / (s^3+1585s^2+1.032e04s+2.351e04
% sysPID = Transfer Function of the MLS System with the PID
sysPID = feedback(mult*kd*sysu*pd,h)

% PID Transfer Function
PIDtf = (kdmult*pd*pic)

% step(sysPID) = Step Response of the Compensated System with the PID
% Peak amplitude = 1.13, Overshoot (%) = 12.8, At time (seconds) = 0.27
% Settling time (seconds) = 0.996
% Rise time (seconds) = 0.00136
step(sysPID),title('Step Response of the Compensated System with the PID')

% General PID Controller Transfer Function
% = (R4/R3)*(R2/R1)*(((R1*C1*s+1)*(R2*C2*s+1))/(R2*C2*s))

% Cancel both R1s and factor out R2 and C2 by multiplying by R2*C2
% = (R4/R3)*R2*(((C1*s+1)*(R2*C2*s+1))/s)

% Overall gain from General PID Controller Transfer Function Equation
% K = (R4*R2*C1)/R3
% Upon the condition that: (R4*R2*C1)/R3 == 167.4176 from kdmult or total
% overall gain of the PID controller.

% PI from General PID Controller Transfer Function Equation
% PI = (s+(1/(R1*C1)))/s
% Upon the condition that: (1/(R1*C1) == 4 from pic transfer function

% PD from General PID Controller Transfer Function Equation
% PD = s+(1/(R2*C2))
% Upon the condition that: (1/(R2*C2) == 3.708 from pd transfer function

% C1 = C2 = 10uF = Engineering Judgement = Often is the case that setting
% capacitor one equal to capacitor two in the General PID Controller
% Transfer Function circuit configuration creates a better scenario for
% involving less variables into the design of the controller.

% From the PI part of the General PID Controller Transfer Function
% Equation:
% (1/(R1*C1) == 4 is equivalent to R1 == 1/(4*C1)

% R1 = 1/(4*10uF) = 1/(4*0.00001) = 25000 = 25 Kilo Ohms

% From the PD part of the General PID Controller Transfer Function
% Equation:
% (1/(R2*C2) == 3.708 is equivalent to R2 == 1/(3.708*C2)

% R2 = 1/(3.708*10uF) = 1/(3.708*0.00001) = 26969 = 26.97 Kilo Ohms

% From the K, or gain, part of the General PID Controller Transfer Function
% Equation:
% (R4*R2*C1)/R3, When R2 = R3 = 26.97 kilo ohms, R2 and R3 cancel and,
% (R4*C1) == 167.4176 is equivalent to R4 == (167.4176/C1)

% R3 = 26.97 Kilo Ohms
% R4 = (167.4176/10uF) = (167.4176/0.00001) = 16741760 = 16.74 Mega Ohms