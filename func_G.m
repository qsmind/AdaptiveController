function b = func_G( sigma_e )
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I_3=eye(3);
b=(1/4)*((1-sigma_e'*sigma_e)*I_3+2*func_S(sigma_e)+2*sigma_e*sigma_e');

end

