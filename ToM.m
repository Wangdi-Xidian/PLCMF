function [C] = ToM(A,row,clom)
%TOM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
C=zeros(row,clom);
for i=1:clom
   C(A(i),i)=1;
end
end

