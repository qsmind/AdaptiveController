% N=10;
% m=1:2:10;
% for i=m
%     f(i)=i;
% end
% f(m)
m=eye(2);
list=1:2:10;
f=[];
for i=list
    f=[f m*i]
end
% f(list);
% i=1:10:N
