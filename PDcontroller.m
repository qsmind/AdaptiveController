steptime=0.02; % 单步时间
start_t=0;  %开始时间
end_t=200;   %结束时间
 %保存结果
result_phi_e_1=[]; 
result_phi_e_2=[]; 
result_phi_e_3=[]; 

result_omega_e_1=[];
result_omega_e_2=[];
result_omega_e_3=[];

result_tau_1=[];
result_tau_2=[];
result_tau_3=[];

result_r_e_1=[];
result_r_e_2=[];
result_r_e_3=[];

result_v_e_1=[];
result_v_e_2=[];
result_v_e_3=[];

result_f_1=[];
result_f_2=[];
result_f_3=[];

result_d_m_tip=[];

result_alpha_tip=[];


result_theta_tri_tip_1=[];
result_theta_tri_tip_2=[];
result_theta_tri_tip_3=[];
result_theta_tri_tip_4=[];
result_theta_tri_tip_5=[];
result_theta_tri_tip_6=[];
result_theta_tri_tip_7=[];

result_rho_tip_1=[];
result_rho_tip_2=[];
result_rho_tip_3=[];


%定义初始值
% t=start_t;
rho=[0.03,0.02,0.025]'; 
m=58.2;

J=[598.3 -22.5 -51.5; 
    -22.5 424.4 -27; 
    -51.5 -27 263.6;];
J_t=[3336.3 -135.4 -154.2;
    -135.4 3184.5 -148.5;
    -154.2 -148.5 2423.7;];
% d_tau=[1+sin(pi*t/125)+sin(pi*t/200) 1+sin(pi*t/125)+sin(pi*t/250) 1+cos(pi*t/125)+cos(pi*t/250)]'*1e-5;
% d_f=[1+sin(pi*t/125)+sin(pi*t/200) 1+sin(pi*t/125)+sin(pi*t/250) 1+cos(pi*t/125)+cos(pi*t/250)]'*1e-4;

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
Lambda_1=0.05*I_3;
Lambda_2=0.2*I_3;

k=0.1;
gamma_1=0.02;
gamma_2=0.02;
gamma_3=0.02;
gamma_4=0.02;
theta_tri_tip=zeros(7,1);
rho_tip=zeros(3,1);
alpha_tip=0;
d_m_tip=0;

v_1=[0.5 0.5 0.5 15 15 15];
K_p=diag(v_1);
v_2=[2 2 2 200 200 200];
K_d=diag(v_2);
%===============

%控制器u
I_6=eye(6);
O_3=zeros(3,3);
O_3x6=zeros(3,6);
O_3x1=zeros(3,1);

% R=I_3-4*(1-sigma_e'*sigma_e)*func_S(sigma_e)/(1+sigma_e'*sigma_e)^2 ...
%     +8*func_S(sigma_e)*func_S(sigma_e)/(1+sigma_e'*sigma_e)^2;
% 
% g=func_S(omega)*v_e+(func_S(omega-omega_e))^2*R*p_t;
% 
% Y_1=g+diag_1*v_e+diag_1*func_S(omega)*r_e;
% Y_2=func_L(func_S(omega_e)*(omega-omega_e))- ...
%     func_S(omega)*func_L(omega)-func_L(diag_2*func_G(sigma_e)*omega_e);
% Y=[Y_1 O_3x6;O_3x1 Y_2];
% 
% C_1=[-func_S(omega) O_3;O_3 O_3];
% C_2=[I_3 O_3;O_3 func_G(sigma_e)];

e_1=[r_e' sigma_e']';
e_2=[v_e' omega_e']';

% diag=[diag_1 O_3;O_3 diag_2];
% s=e_2+diag*e_1;
% A_tip=[O_3 O_3;func_S(rho_tip) O_3;];
% K=[K_1 O_3;O_3 K_2];

%==================================

%对象

% I_3=eye(3);

% sigma_e=(sigma_t*(sigma'*sigma-1)+sigma*(1-sigma_t'*sigma_t)-2*func_S(sigma_t)*sigma)...
%     /(1+sigma_t'*sigma_t*sigma'*sigma+2*sigma_t'*sigma);

% d=[d_f;d_tau];
% M=[m*I_3 O_3;O_3 J];
% h=[-m*g;-func_S(omega)*J*omega-J*func_S(omega)*omega_e];
% N=[-m*I_3 O_3;O_3 J];
% B=[R*func_S(p_t);R];
% A=[O_3 O_3;func_S(rho) O_3];
% eta=R'*(omega-omega_e);

% ====================

% 迭代
% start_t=start_t+steptime;
for t=start_t:steptime:end_t

d_tau=[1+sin(pi*t/125)+sin(pi*t/200) ;
    1+sin(pi*t/125)+sin(pi*t/250) ;
    1+cos(pi*t/125)+cos(pi*t/250)]*1e-5;
d_f=[1+sin(pi*t/125)+sin(pi*t/200);
    1+sin(pi*t/125)+sin(pi*t/250);
    1+cos(pi*t/125)+cos(pi*t/250)]*1e-4;


% 控制器
R=I_3-4*(1-sigma_e'*sigma_e)*func_S(sigma_e)/(1+sigma_e'*sigma_e)^2 ...
    +8*func_S(sigma_e)*func_S(sigma_e)/(1+sigma_e'*sigma_e)^2;
% R=I_3-(4*(1-sigma_e'*sigma_e)/(1+sigma_e'*sigma_e)^2)*func_S(sigma_e) ...
%     +(8/(1+sigma_e'*sigma_e)^2)*func_S(sigma_e)*func_S(sigma_e);

g=func_S(omega)*v_e+(func_S(omega-omega_e))^2*R*p_t;

Y_1=g+Lambda_1*v_e+Lambda_1*func_S(omega)*r_e;
Y_2=func_L(func_S(omega_e)*(omega-omega_e))- ...
    func_S(omega)*func_L(omega)-func_L(Lambda_2*func_G(sigma_e)*omega_e);
Y=[Y_1 O_3x6;
    O_3x1 Y_2];

C_1=[-func_S(omega) O_3;
    O_3 O_3];
C_2=[I_3 O_3;
    O_3 func_G(sigma_e)];

Lambda=[Lambda_1 O_3;
    O_3 Lambda_2];
s=e_2+Lambda*e_1;
A_tip=[O_3 O_3;
    func_S(rho_tip) O_3;];
K=[K_1 O_3;
    O_3 K_2];

u=-K_p*e_1-K_d*e_2;

f=u(1:3,1);
tau=u(4:6,1);
%==================================

%对象

d=[d_f;
    d_tau];
M=[m*I_3 O_3;
    O_3 J];
h=[-m*g;
    -func_S(omega)*J*omega-J*func_S(omega)*omega_e];
N=[-m*I_3 O_3;
    O_3 J];
B=[R*func_S(p_t);
    R];
A=[O_3 O_3;
    func_S(rho) O_3];

eta=R'*(omega-omega_e);

de_1=C_1*e_1+C_2*e_2;
de_2=(M^-1)*(h+N*B*(J_t^-1)*func_S(eta)*J_t*eta+(I_6+A)*u+d);

% =================================================


%参数变化率
dtheta_tri_tip=gamma_1*Y'*s;
drho_tip=-gamma_2*(func_S(f))'*(omega_e+Lambda_2*sigma_e);
dalpha_tip=gamma_3*(1+norm(p_t,2))*norm(omega-omega_e,2)^2*norm(s,1);
dd_m_tip=gamma_4*norm(s,1);


% 转为显示角度
phi_e_1=atand(R(1,2)/R(1,1));
phi_e_2=-asind(R(1,3));
phi_e_3=atand(R(2,3)/R(3,3));

% 保存结果
result_phi_e_1=[result_phi_e_1,phi_e_1];
result_phi_e_2=[result_phi_e_2,phi_e_2];
result_phi_e_3=[result_phi_e_3,phi_e_3];

result_omega_e_1=[result_omega_e_1,omega_e(1,1)];
result_omega_e_2=[result_omega_e_2,omega_e(2,1)];
result_omega_e_3=[result_omega_e_3,omega_e(3,1)];

result_tau_1=[result_tau_1,tau(1,1)];
result_tau_2=[result_tau_2,tau(2,1)];
result_tau_3=[result_tau_3,tau(3,1)];

result_r_e_1=[result_r_e_1,r_e(1,1)];
result_r_e_2=[result_r_e_2,r_e(2,1)];
result_r_e_3=[result_r_e_3,r_e(3,1)];

result_v_e_1=[result_v_e_1,v_e(1,1)];
result_v_e_2=[result_v_e_2,v_e(2,1)];
result_v_e_3=[result_v_e_3,v_e(3,1)];

result_f_1=[result_f_1,f(1)];
result_f_2=[result_f_2,f(2)];
result_f_3=[result_f_3,f(3)];

% ===============================================================


% 更新参数
theta_tri_tip=func_update(theta_tri_tip,dtheta_tri_tip,steptime);
rho_tip=func_update(rho_tip,drho_tip,steptime);
alpha_tip=func_update(alpha_tip,dalpha_tip,steptime);
d_m_tip=func_update(d_m_tip,dd_m_tip,steptime);

e_1=func_update(e_1,de_1,steptime);
e_2=func_update(e_2,de_2,steptime);

r_e=e_1(1:3,1);
sigma_e=e_1(4:6,1);
v_e=e_2(1:3,1);
omega_e=e_2(4:6,1);

% ============================
end

% 画图
t=start_t:steptime:end_t;

% ================================================================================================
figure(1)
%绘制phi
% figure
subplot(3,1,1)
plot(t,result_phi_e_1,':',t,result_phi_e_2,'--',t,result_phi_e_3,'-')
% ylim([-50 50])
% xlabel('t(s)')
ylabel('$ \Phi_e $(deg) ','Interpreter','latex','Fontsize',12)
legend({'$ {\Phi_e}_1 $','$ {\Phi_e}_2 $','$ {\Phi_e}_3 $'},'Interpreter','latex','Fontsize',12)

%绘制omega
% figure
subplot(3,1,2)
plot(t,result_omega_e_1,':',t,result_omega_e_2,'--',t,result_omega_e_3,'-')
% ylim([-50 50])
% xlabel('t(s)')
ylabel('$\omega_e $(rad/s)','Interpreter','latex','Fontsize',12)
legend({'${\omega_e}_1$','${\omega_e}_2$','${\omega_e}_3$'},'Interpreter','latex','Fontsize',12)

%绘制tau
% figure
subplot(3,1,3)
plot(t,result_tau_1,':',t,result_tau_2,'--',t,result_tau_3,'-')
% ylim([-50 50])
xlabel('$t(s)$','Interpreter','latex','Fontsize',12)
ylabel('$\tau $(Nm)','Interpreter','latex','Fontsize',12)
legend({'$\tau_1$','$\tau_2$','$\tau_3$'},'Interpreter','latex','Fontsize',12)
% ================================================================================================
figure(2)
%绘制r_e
% figure
subplot(3,1,1)
plot(t,result_r_e_1,':',t,result_r_e_2,'--',t,result_r_e_3,'-')
% ylim([-50 50])
% xlabel('t(s)')
ylabel('$r_e$(m)','Interpreter','latex','Fontsize',12)
legend({'${r_e}_1$','${r_e}_2$','${r_e}_3$'},'Interpreter','latex','Fontsize',12)

%绘制v_e
% figure
subplot(3,1,2)
plot(t,result_v_e_1,':',t,result_v_e_2,'--',t,result_v_e_3,'-')
% ylim([-50 50])
% xlabel('t(s)')
ylabel('$v_e$(m/s)','Interpreter','latex','Fontsize',12)
legend({'${v_e}_1$','${v_e}_2$','${v_e}_3$'},'Interpreter','latex','Fontsize',12)

%绘制f
% figure
subplot(3,1,3)
plot(t,result_f_1,':',t,result_f_2,'--',t,result_f_3,'-')
% ylim([-50 50])
xlabel('$t(s)$','Interpreter','latex','Fontsize',12)
ylabel('$ f $(N) ','Interpreter','latex','Fontsize',12)
legend({'$f_1$','$f_2$','$f_3$'},'Interpreter','latex','Fontsize',12)
% ================================================================================================
