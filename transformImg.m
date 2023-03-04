%
function P = transformImg(img,M)

%i -> column
%j -> row. both starting with zero

P=zeros(size(img,1)*size(img,2),4);
c=0;
for i=1:size(img,2)
    for j=1:size(img,1)
        Pxyz =  double(M) * double([i-1;j-1;0;1]);
        c=c+1;
        P(c,1:4) = [Pxyz(1:3)' double(img(j,i)) ];   %location value
    end
end
