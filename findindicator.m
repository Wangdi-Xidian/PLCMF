function outVec = findindicator(xVec,C)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
c=size(C,2);
obj = zeros(c, 1);
tmp = eye(c);
for i=1:c
    obj(i,1) = obj(i,1) + (norm(xVec - C(:,i))^2);
end
[min_val, min_idx] = min(obj);
outVec = tmp(:, min_idx);
end


