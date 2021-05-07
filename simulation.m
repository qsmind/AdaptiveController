% 单步时间
steptime=0.01;
%定义初始值
t=0;
rho=[0.03,0.02,0.025]'; 
m=58.2;

J=[598.3 -22.5 -51.5;-22.5 424.4 -27;-51.5 -27 263.6;];
J_t=[3336.3 -135.4 -154.2;-135.4 3184.5 -148.5;-154.2 -148.5 2423.7;];
d_tau=[1+sin(pi*t/125)+sin(pi*t/200) 1+sin(pi*t/125)+sin(pi*t/250) 1+cos(pi*t/125)+cos(pi*t/250)]'*1e-5;
d_f=[1+sin(pi*t/125)+sin(pi*t/200) 1+sin(pi*t/125)+sin(pi*t/250) 1+cos(pi*t/125)+cos(pi*t/250)]'*1e-4;

theta_0=[50,580,400,300,0,0,0]';
p_t=[0 5 0]';

r=[1 1 1]'*7.078*1e6;
v=[0 0 0]';
sigma=[0 0 0]';
omega=[0 0 0]';
r_e=[50/sqrt(2) 0 -50/sqrt(2)]';
v_e=[0.5 -0.5 0.5]';
sigma_e=[0.5 -0.6 0.7]';
omega_e=[0.02 -0.02 0.02]';

I_3=eye(3);
K_1=20*I_3;
K_2=200*I_3;
diag_1=0.05*I_3;
diag_2=0.2*I_3;
k=0.1;
gamma_1=0.02;
gamma_2=0.02;
gamma_3=0.02;
gamma_4=0.02;
theta_tri_tip=zeros(7,1);
rho_tip=zeros(3,1);
alpha_tip=0;
d_m_tip=0;
%===============



%控制器u
I_6=eye(6);
O_3=zeros(3,3);
O_3x6=zeros(3,6);
O_3x1=zeros(3,1);

R=I_3-4*(1-sigma_e'*sigma_e)*func_S(sigma_e)/(1+sigma_e'*sigma_e)^2 ...
    +8*func_S(sigma_e)*func_S(sigma_e)/(1+sigma_e'*sigma_e)^2;

g=func_S(omega)*v_e+(func_S(omega-omega_e))^2*R*p_t;

Y_1=g+diag_1*v_e+diag_1*func_S(omega)*r_e;
Y_2=func_L(func_S(omega_e)*(omega-omega_e))- ...
    func_S(omega)*func_L(omega)-func_L(diag_2*func_G(sigma_e)*omega_e);
Y=[Y_1 O_3x6;O_3x1 Y_2];

C_1=[-func_S(omega) O_3;O_3 O_3];
C_2=[I_3 O_3;O_3 func_G(sigma_e)];

e_1=[r_e' sigma_e']';
e_2=[v_e' omega_e']';

r_e=e_1(1:3,1);
sigma_e=e_1(4:6,1);
v_e=e_2(1:3,1);
omega_e=e_2(4:6,1);

diag=[diag_1 O_3;O_3 diag_2];
s=e_2+diag*e_1;
A_tip=[O_3 O_3;func_S(rho_tip) O_3;];
K=[K_1 O_3;O_3 K_2];
u=((I_6+A_tip)^-1)*(-k*C_2'*e_1-K*s-Y*(theta_0+theta_tri_tip)- ...
    d_m_tip*sign(s)-alpha_tip*(1+norm(p_t,2))*norm(omega-omega_e)^2*sign(s));


f=u(1:3,1);
tau=u(4:6,1);
%==================================

%对象

I_3=eye(3);

% sigma_e=(sigma_t*(sigma'*sigma-1)+sigma*(1-sigma_t'*sigma_t)-2*func_S(sigma_t)*sigma)...
%     /(1+sigma_t'*sigma_t*sigma'*sigma+2*sigma_t'*sigma);


d=[d_f;d_tau];
M=[m*I_3 O_3;O_3 J];
h=[-m*g;-func_S(omega)*J*omega-J*func_S(omega)*omega_e];
N=[-m*I_3 O_3;O_3 J];
B=[R*func_S(p_t);R];
A=[O_3 O_3;func_S(rho) O_3];

de_1=C_1*e_1+C_2*e_2;
de_2=(M^-1)*(h+N*B*(J_t^-1)*func_S(eta)*J_t*eta+(I_6+A)*u+d);

eta=R'*(omega-omega_e);

% ====================

%参数变化率
dtheta_tri_tip=gamma_1*Y'*s;
drho_tip=-gamma_2*(func_S(f))'*(omega_e+diag_2*sigma_e);
dalpha_tip=gamma_3*(1+norm(p_t,2))*norm(omega-omega_e)^2*norm(s,1);
dd_m_tip=gamma_4*norm(s,1);

% 更新参数
theta_tri_tip=func_update(theta_tri_tip,dtheta_tri_tip,steptime);
rho_tip=func_update(rho_tip,drho_tip,steptime);
alpha_tip=func_update(alpha_tip,dalpha_tip,steptime);
d_m_tip=func_update(d_m_tip,dd_m_tip,steptime);

e_1=func_update(e_1,de_1,steptime);
e_2=func_update(e_2,de_2,steptime);

% ============================
phi_e_1=atan(R(1,2)/R(1,1));
phi_e_2=-asin(R(1,3));
phi_e_3=atan(R(2,3)/R(3,3));

