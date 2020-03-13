% code to separate the rays into groups


for i=1:60000*20
    r(i).l=i;
    r(i).j=i+67;
end
tic()
b=[1,2,3,4]; %vettore dei numeri
vec= {r.l}';

vec1=cell2mat(vec);
sizeb=4; %ho 4 gruppi
for i=1:sizeb
    A=find(vec1(:,1)==b(i));
    eval(['Ray_' num2str(i) '=r(A);']); %smisto
    
end


toc()
A=1;
% a=1; b=2; c=3;
% r(300000).m=0;
% v=zeros(1,20);
% v(1,1)=1;
% v(1,2)=1;
% 
% tic()
% for i=1:300000
%    A=find(v(1,:)==0);
%    p=v;
%    p(1,A(1))=2;
%    r(i).m=p;
%    r(i).ap=33;
%    %r(i).m=sscanf(sprintf('%d',[a b c]),'%d');
%     % r(i).m=str2num(sprintf('%-1d',[a,b,c]));
% end
% toc()

toc()

%
%
a(1).m=[1,2,3,4,5];
a(2).m=[1,2,3,0,0];
a(3).m=[1,2,3,0,0];
a(4).m=[1,2,3,4,5];
a(5).m=[1,2,0,0,0];
a(6).m=[1,0,0,0,0];
memory={a.m}';
r(1).m=1;
r(2).m=2;
r(3).m=3;
r(4).m=4;
r(5).m=5;
r(6).m=6;
memo=cell2mat(memory);
Groups=unique(memo, 'rows');
[g,j]=size(Groups);

for i=1:g
  %  A=find(memo(:,:)==me(i,:));
    tf = find((ismember(memo, Groups(i,:), 'rows')));
    eval(['Ray_' num2str(i) '=r(tf);']); %smisto 
end


a=1; b=2; c=3;
r(300000).m=0;
v=zeros(1,20);
v(1,1)=1;
v(1,2)=1;
tic()
for i=1:300000
   idx = find(v(1,:)==0, 1, 'first');
  % A=find(v(1,:)==0);
  % p=v;
%    p(1,A(1))=2;
%    r(i).m=p;
   r(i).ap=33;
   % r(i).m=sscanf(sprintf('%d',[a b c]),'%d');
    % r(i).m=str2num(sprintf('%-1d',[a,b,c]));
end
toc()






