% X=imread('2021-0705-number3 sunfengqin 64 years old-Actin&Perforin_Series001_Lng_Processed001.jpg');
% [m n p]=size(X);
% for i=1:m
%     for j=1:n
%         X_red(:,:)=X(:,:,1);
%         X_green(:,:)=X(:,:,2);
%         X_blue(:,:)=X(:,:,3);
%     end
% end

Y=imread('VS1.tif');
[m n]=size(Y);
distance=zeros(50000,1);
p=zeros(50000,1);
Xaxis=zeros(50000,1);
Yaxis=zeros(50000,1);
Gray=zeros(50000,1);
a=0;%a用于计算高于阈值荧光像素点计数
x1=input('请输入x1=');
y1=input('请输入y1=');
x2=input('请输入x2=');
y2=input('请输入y2=');
k=(y2-y1)/(x2-x1);
b=y1-k*x1;
x3=input('请输入x3=');
y3=input('请输入y3=');
z=(abs(k*x3-y3+b)/sqrt(k*k+1));
for i=1:m
    for j=1:n
        if Y(i,j)>10
        a=a+1;
        Xaxis(a)=j;
        Yaxis(a)=i;
        distance(a)=(abs(k*j-i+b)/sqrt(k*k+1));
        Gray(a)=Y(i,j);
        if distance(a)<z
        p(a)=distance(a)/z;
        end
        end
        
        
    end
end
s=sum(Gray);
E=zeros(10,1);
for t=1:a
    if p(t)>0
    g=floor(p(t)/0.2)+1;
    E(g,1)=Gray(t,1)+E(g,1);
    end
end
E=E./s;
C=zeros(10,1);
d=sum(E);
C=E./d;
