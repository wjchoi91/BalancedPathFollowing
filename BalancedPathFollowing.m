clear all

%% Pricipal particulars of a KVLCC2 tanker
Lpp = 320.00;         B = 58;           d = 20.8;
disp = 312600;        xG = 11.2;         Cb = 0.810;
DP = 9.86;        HR = 15.8;        AR = 112.5;

%% Hydrodynamic force coefficient
R0_ndim = 0.022;       Xvv_ndim = -0.040;     Xvr_ndim = 0.002;
Xrr_ndim = 0.011;      Xvvvv_ndim = 0.771;    Yv_ndim = -0.315;
YR_ndim = 0.083;       Yvvv_ndim = -1.607;    Yvvr_ndim = 0.379;
Yvrr_ndim = -0.391;    Yrrr_ndim = 0.008;     Nv_ndim = -0.137;
NR_ndim = -0.049;      Nvvv_ndim = -0.030;    Nvvr_ndim = -0.294;
Nvrr_ndim = 0.055;     Nrrr_ndim = -0.013;
mx_ndim = 0.022;       my_ndim = 0.223;       Jz_ndim = 0.011;
tP = 0.220;        tR = 0.387;        aH = 0.312;
xH_ndim = -0.464;      C1 = 2.0;
C21 = 1.6;                 % beta_P > 0
C22 = 1.1;                 % beta_P < 0
gamma_R1 = 0.640;           % beta_R > 0
gamma_R2 = 0.395;           % beta_R < 0
lR_ndim = -0.710;      epsilon = 1.09;     kappa = 0.50;
f_alpha = 2.747;


%% Other coefficient
k0 = 0.293;       k1 = -0.275;      k2 = -0.139;
wP0 = 0.35;        nP = 1.53;        xP_ndim = -0.48;
rho = 1025;         eta = DP / HR;    xR = -0.5*Lpp;
mass_ndim_coef = 0.5*rho*Lpp^2*d;
moment_ndim_coef = 0.5*rho*Lpp^4*d;
xH = xH_ndim*Lpp;
m = disp*rho;
mx = mx_ndim * mass_ndim_coef;
my = my_ndim * mass_ndim_coef;
Jz = Jz_ndim * moment_ndim_coef;
IzG = (0.25*Lpp)^2*disp*rho;

%% Initialize
u = 7.9732;         vm = 0;             r = 0;
Accel = [0; 0; 0]; % u_dot, vm_dot, r_dot
Velocity = [u; 0; 0]; % u, vm, r
Position = [0; 0; 0]; % x0, y0, psi
u_dot = 0;          vm_dot = 0;         r_dot = 0;
delta = 0;
KP = 8.2;             KD = 173.2;
psi_order = 0;
Time_s = 1;
Time_f = 2093;

Advance35_Stbd = 1038.2;
Transfer35_Stbd = 438.5;
Advance35_Port = 988.9;
Transfer35_Port = 396.4;

%% Waypoint
WP_List = [0 0;
    2600 0;
    4765 1250;
    7663 473;
    8388 -2232;
    5973 -2879];
WP = [WP_List(2,1), WP_List(2,2)];

Position(1) = 200;
Position(2) = 0;
Position(3) = pi/2;
curr_WP_No = 2;
WP_change = 1;


for i = 1:(Time_f/Time_s)+1
        
    %% Time
    t = i*Time_s;
        
    %% Position, velocity, angular velocity
    u = Velocity(1);        vm = Velocity(2);       r = Velocity(3);
    u_dot = Accel(1);       vm_dot = Accel(2);      r_dot = Accel(3);
    x0 = Position(1);       y0 = Position(2);       psi = Position(3); 
    U = sqrt(u^2 + vm^2);
    v = vm + xG*r;          vm_ndim = vm/U;         r_ndim = r*Lpp/U;
    beta = atan2(-vm,u);
        
    %% Hydrodynamic forces acting on ship hull
    XH_ndim = -R0_ndim + Xvv_ndim*vm_ndim^2 + Xvr_ndim*vm_ndim*r_ndim + Xrr_ndim*r_ndim^2 + Xvvvv_ndim*vm_ndim^4;
    YH_ndim = Yv_ndim*vm_ndim + YR_ndim*r_ndim + Yvvv_ndim*vm_ndim^3 + Yvvr_ndim*vm_ndim^2*r_ndim + Yvrr_ndim*vm_ndim*r_ndim^2 + Yrrr_ndim*r_ndim^3;
    NH_ndim = Nv_ndim*vm_ndim + NR_ndim*r_ndim + Nvvv_ndim*vm_ndim^3 + Nvvr_ndim*vm_ndim^2*r_ndim + Nvrr_ndim*vm_ndim*r_ndim^2 + Nrrr_ndim*r_ndim^3;
    
    XH = 0.5*rho*Lpp*d*U^2*XH_ndim;
    YH = 0.5*rho*Lpp*d*U^2*YH_ndim;
    NH = 0.5*rho*Lpp^2*d*U^2*NH_ndim;
    
    %% Hydrodynamic force due to propeller
    beta_P = beta - xP_ndim*r_ndim;
    if beta_P >= 0
        C2 = C21;
    else
        C2 = C22;
    end
    wP = wP0*exp(-4*beta_P^2);
    JP = u*(1 - wP) / (nP*DP);
    KT = k2*JP^2 + k1*JP + k0;
    T = rho*nP^2*DP^4*KT;
    
    XP = (1 - tP)*T;
    
    %% Hydrodynamic force by steering
    beta_R = beta - lR_ndim*r_ndim;
    if beta_R >= 0
        gamma_R = gamma_R1;
    else
        gamma_R = gamma_R2;
    end
    vR = U*gamma_R*beta_R;
    uR = epsilon*u*(1 - wP) * sqrt(eta*(1 + kappa*(sqrt(1 + (8*KT)/(pi*JP^2)) - 1))^2 + (1 - eta));
    alpha_R = delta - vR / uR;
    UR = sqrt(uR^2 + vR^2);
    FN = 0.5*rho*AR*UR^2*f_alpha*sin(alpha_R);
    
    XR = -(1 - tR)*FN*sin(delta);
    YR = -(1 + aH)*FN*cos(delta);
    NR = -(xR + aH*xH)*FN*cos(delta);
    
    %% Motion equation
    X = XH + XR + XP;
    Y = YH + YR;
    Nm = NH + NR;

    M = [m+mx 0 0; 0 m+my xG*m; 0 xG*m IzG+xG^2*m+Jz];
    C = [0 -(m+my)*r -xG*m*r; (m+mx)*r 0 0; xG*m*r 0 0];
    P = [X; Y; Nm];
    Vel = Velocity;

    Accel = M\(P-C*Vel);
    
    Velocity = Velocity + Accel*Time_s;
    J = [cos(psi) -sin(psi) 0; sin(psi) cos(psi) 0; 0 0 1];
    Position_dot = J * Velocity;
    Position = Position + Position_dot*Time_s;
    
    %% Waypoint swtiching
    if curr_WP_No ~= size(WP_List, 1)
        if WP_change == 1
            Curr_Co = atan2(WP_List(curr_WP_No,1)-WP_List(curr_WP_No-1,1), WP_List(curr_WP_No,2)-WP_List(curr_WP_No-1,2));
            Next_Co = atan2(WP_List(curr_WP_No+1,1)-WP_List(curr_WP_No,1), WP_List(curr_WP_No+1,2)-WP_List(curr_WP_No,2));
            Alt_Co = (Next_Co - Curr_Co);
            
            if Alt_Co < -pi
                Alt_Co = Alt_Co + 2 * pi;
            elseif Alt_Co > pi
                Alt_Co = Alt_Co - 2 * pi;
            end
            
            if Alt_Co >= 0
                Advance = Advance35_Stbd;
                Transfer = Transfer35_Stbd;
            else
                Advance = Advance35_Port;
                Transfer = Transfer35_Port;
            end
            
            Dist_New_Co = (Advance - ( (Transfer - Transfer * tan((pi/2-abs(Alt_Co))/2)) / tan(abs(Alt_Co)) ))*1;
            WP_change = 0;
        end
        
        WP_dist = sqrt((WP_List(curr_WP_No, 2) - Position(1))^2 + (WP_List(curr_WP_No, 1) - Position(2))^2);
        if WP_dist < Dist_New_Co
            curr_WP_No = curr_WP_No + 1;
            WP_change = 1;
        end
    end

    %% Balanced LOS guidance
    Curr_Co = atan2(WP_List(curr_WP_No,1)-WP_List(curr_WP_No-1,1), WP_List(curr_WP_No,2)-WP_List(curr_WP_No-1,2));
    XTE = -(Position(1)-WP_List(curr_WP_No, 2))*sin(Curr_Co)+(Position(2)-WP_List(curr_WP_No, 1))*cos(Curr_Co);

    Diff_Co = Curr_Co - psi;
        
    if Diff_Co < -pi
        Diff_Co = Diff_Co + 2 * pi;
    elseif Diff_Co > pi
        Diff_Co = Diff_Co - 2 * pi;
    end

    Dist_New_Co2 = (Advance - ( (Transfer - Transfer * tan((pi/2-abs(Diff_Co))/2)) / tan(abs(Diff_Co)) ))*1;

    n1 = Dist_New_Co2+0.25*Lpp^2/abs(XTE);
    if n1 < abs(XTE)
        n1 = abs(XTE);
    end

    % Calculate LOS position
    delta_Lat1 = WP_List(curr_WP_No,2) - WP_List(curr_WP_No-1,2); % x
    delta_Long1 = WP_List(curr_WP_No,1) - WP_List(curr_WP_No-1,1); % y

    if delta_Long1 ~= 0        
        if delta_Lat1 == 0
            LOS_Lat1 = WP_List(curr_WP_No,2);
            if delta_Long1 > 0
                LOS_Long1 = Position(2) + sqrt( (n1)^2 + (WP_List(curr_WP_No,2) - Position(1))^2 );
            elseif delta_Long1 < 0
                LOS_Long1 = Position(2) - sqrt( (n1)^2 + (WP_List(curr_WP_No,2) - Position(1))^2 );
            end
        else
            coef_d1 = delta_Long1 / delta_Lat1;
            coef_e1 = WP_List(curr_WP_No-1,2);
            coef_f1 = WP_List(curr_WP_No-1,1);
            coef_g1 = coef_f1 - coef_d1 * coef_e1;
            coef_a1 = 1 + coef_d1^2;
            coef_b1 = 2 * (coef_d1 * coef_g1 - coef_d1 * Position(2) - Position(1));
            coef_c1 = Position(1)^2 + Position(2)^2 + coef_g1^2 - (n1)^2 - 2 * coef_g1 * Position(2);
            if delta_Lat1 > 0
                LOS_Lat1 = (-coef_b1 + sqrt(coef_b1^2 - 4*coef_a1*coef_c1)) / (2*coef_a1);
            elseif delta_Lat1 < 0
                LOS_Lat1 = (-coef_b1 - sqrt(coef_b1^2 - 4*coef_a1*coef_c1)) / (2*coef_a1);
            end
            LOS_Long1 = coef_d1 * (LOS_Lat1 - coef_e1) + coef_f1;
        end
    else
        LOS_Long1 = WP_List(curr_WP_No, 1);
    end

    % Calculate psi_order
    psi_order = atan2(LOS_Long1 - Position(2), LOS_Lat1 - Position(1));
     
    %% Heaing controller
    error = psi_order - psi;
    if error > pi
        error = error - 2*pi;
    elseif error < -pi
        error = error + 2*pi;
    end

    delta_order = KP*error - KD*r*Time_s;

    if abs(delta_order - delta) <= 2.32*pi/180*Time_s
        delta = delta_order;
    elseif delta_order - delta > 2.32*pi/180*Time_s
        delta = delta + 2.32*pi/180*Time_s;
    elseif delta_order - delta < -2.32*pi/180*Time_s
        delta = delta - 2.32*pi/180*Time_s;
    end
    
    if delta > 35*pi/180
        delta = 35*pi/180;
    elseif delta < -35*pi/180
        delta = -35*pi/180;
    end
    
    %% Store data
    Data(i+1, 1) = t;
    Data(i+1, 2:4) = Position';
    Data(i+1, 5:7) = Velocity';
    Data(i+2, 8) = delta;
    Data(i+1, 9) = U;
    Data(i+1, 10) = n1;

end

figure(1)
subplot(3, 2, [1,3,5]);
plot(Data(2:end-1,3), Data(2:end-1,2), 'r', 'LineWidth', 3);
hold on
plot(WP_List(:,1), WP_List(:,2), '--ob', 'LineWidth', 3)
hold off
axis('equal');
grid on;
axis('equal');
xlabel('y_0 [m]');
ylabel('x_0 [m]');
subplot(3, 2, 2)
plot(Data(2:end-1,1), Data(2:end-1,9), 'b', 'LineWidth', 3);
grid on;
xlabel('Time [sec]');
ylabel('u [m/s]');
subplot(3, 2, 4)
plot(Data(2:end-1,1), Data(2:end-1,8)*180/pi, 'b', 'LineWidth', 3);
grid on;
xlabel('Time [sec]');
ylabel('delta [deg]');
subplot(3, 2, 6)
plot(Data(2:end-1,1), Data(2:end-1,7)*180/pi, 'b', 'LineWidth', 3);
grid on;
xlabel('Time [sec]');
ylabel('r [deg/s]');