function   [result]=PLCMF(X,gt,option)
%%----------------Initialize-------------------
numClust=option.numClust ;
K=option.K;
threshold=option.threshold;
delta=option.delta;
lambda=option.lambda;         
r=option.r;
max_iter=option.max_iter;
Vnum=option.Vnum;
N=option.N ;                            
alpha=option.alpha; 
beta=option.beta; 
alpha_r=alpha.^r;
Jlast=99999;
IsConverge = 0;
iter = 1;
Q=cell(size(X,1),size(X,2));
U=cell(size(X,1),size(X,2));
Y=cell(size(X,1),size(X,2));
V1=zeros(K,K);
V2=zeros(K,N);
V=rand(K,N);
J1=ones(Vnum,1);
J2=ones(Vnum,1);
J4=ones(Vnum,1);
Ja=ones(Vnum,1);
tic;
for i=1:Vnum 
    Q{i}=rand(numClust,K);
    U{i}=rand(size(X{i},1),K);
    Y{i}=rand(numClust,N);
    YY=litekmeans(X{i}',numClust,'MaxIter', 100);
    Y{i}=ToM(YY,size(Y{i},1),size(Y{i},2));
end 
C=rand(K,numClust);
g0=ceil(rand(1,N)*numClust);
G=ToM(g0',numClust,size(V,2));                          
%%------------------Update---------------------------

while (IsConverge == 0&&iter<max_iter+1) 
    %---------UpdateQi---------
     for i=1:Vnum
        Q{i}=Y{i}*V'/(V*V'+lambda/(alpha_r(i)*delta)*eye(size(V,1)));
     end
     %---------UpdateUi---------
     for i=1:Vnum
        U{i}=X{i}*V'/(V*V'+lambda/(alpha_r(i))*eye(size(V,1)));
     end     
    %---------UpdateV---------
    V1=zeros(K,K);
    V2=zeros(K,N);
     for i=1:Vnum
        V1=V1+(alpha_r(i)*U{i}'*U{i}+alpha_r(i)*delta*Q{i}'*Q{i});
        V2=V2+(alpha_r(i)*U{i}'*X{i}+alpha_r(i)*delta*Q{i}'*Y{i});
     end  
     V=(V1+(delta+lambda)*eye(K))\(V2+beta*C*G);   

     %---------UpdateC---------
     C=V*G'/(G*G'+lambda/beta*eye(size(numClust)));
     %---------UpdateG---------
     for i = 1:N
        xVec = V(:,i);
        G(:,i) = findindicator(xVec, C);
     end

    %---------Calculate J---------
    J3=norm(V-C*G,'fro')^2;
    for i=1:Vnum 
        J1(i)=alpha_r(i)*norm(X{i}-U{i}*V,'fro')^2;
        J2(i)=alpha_r(i)*norm(Y{i}-Q{i}*V,'fro')^2;
        J4(i)=norm(U{i},'fro')^2+norm(Q{i},'fro')^2;
        Ja(i)=J1(i)+delta*J2(i);   
    end
    J5=(norm(V,'fro')^2+norm(C,'fro')^2)+sum(J4(i));
    Jcurrent=sum(J1)+delta*sum(J2)+beta*J3+lambda*J5;
    %---------Calculate alpha---------   
    alpha=(Ja.^(1/(1-r)))/sum(Ja.^(1/(1-r)));
    alpha_r=alpha.^r;
    
    %---------Iscoverage---------   
    J(iter)=Jcurrent;
    if (abs(Jlast - Jcurrent)) < threshold
            IsConverge=1;
    end
    Jlast=Jcurrent;
    iter=iter+1; 
end
toc
[~,pred_label] = max(G,[],1);   
pred_label=pred_label';
 result = ClusteringMeasure(gt, pred_label);                           
 [f,p,r] = compute_f(gt,pred_label);                                   
 [A nmi avgent] = compute_nmi(gt,pred_label);
 result=[result f p r];
end


