function [C] = ToM(A,row,clom)
%TOM 此处显示有关此函数的摘要
%   此处显示详细说明
C=zeros(row,clom);
for i=1:clom
   C(A(i),i)=1;
end
end

