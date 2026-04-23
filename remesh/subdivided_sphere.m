

function[VV,FF] = subdivided_sphere(N)
%Icosahedron coordinates
%-----------------------------------
t=(1+sqrt(5))/2;
V=1/sqrt(1+t^2)*[t 1 0; -t 1 0; t -1 0; -t -1 0; 1 0 t; 1 0 -t; -1 0 t; -1 0 -t; 0 t 1; 0 -t 1; 0 t -1; 0 -t -1];
F=[1 9 5; 1 6 11; 3 5 10; 3 12 6; 2 7 9; 2 11 8; 4 10 7; 4 8 12; 1 11 9; 2 9 11; 3 10 12; 4 12 10; 5 3 1; 6 1 3; 7 2 4; 8 4 2; 9 7 5; 10 5 7; 11 6 8; 12 8 6];

VV=[];
FF=[];

%Divide each face into triangles
%-----------------------------------
for i=1:1:size(F,1)
    
v1=V(F(i,1),:);
v2=V(F(i,2),:);
v3=V(F(i,3),:);

%Vertices of the triangles
c=0;
for j=0:1:N
v12=(v1*(N-j)+v2*j)/N;
v13=(v1*(N-j)+v3*j)/N;
for k=0:1:j
c=c+1;    
if j==0
v123(c,:)=v12;
else
v123(c,:)=(v12*(j-k)+v13*k)/j;    
end
    
end
         
end
 
%Faces of the triangles
c=0;
FT=[];

for j=1:1:N
for k=1:1:j
c=c+1;    
n=j*(j-1)/2+k;
FT(c,:)=[n n+j n+j+1];    
end  
end

for j=2:1:N
for k=1:1:j-1
c=c+1;    
n=j*(j-1)/2+k;
FT(c,:)=[n n+j+1 n+1];    
end  
end

%Merge vertices and faces
m=size(VV,1);
VV=[VV; v123];
FF=[FF; FT+m];
       
end

%Remove duplicate vertices
%-----------------------------------
q=size(VV,1);
q2=q;
for i=1:1:q-1
for j=i+1:q    
if norm(VV(i,:)-VV(j,:))<1e-5
VV(j,:)=[];
q2=q2-1;
for k=1:1:size(FF,1)
for l=1:1:size(FF,2)
if FF(k,l)==j 
FF(k,l)=i;
elseif FF(k,l)>j
FF(k,l)=FF(k,l)-1;    
end    
end
end    
end
if j>=q2; break; end
end
if i>=q2-1; break; end
end

%Project vertices into unit sphere
%-----------------------------------
for i=1:1:q2
VV(i,:)=VV(i,:)/norm(VV(i,:));       
end

%Plot
%-----------------------------------

end