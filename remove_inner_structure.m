function [p,eee,tt]=remove_inner_structure(p,e,t,n)
p1=p(:,t(1,:));
p2=p(:,t(2,:));
p3=p(:,t(3,:));
m=(p1+p2+p3)/3; % midpoints
mx=floor(n*m(1,:));
my=floor(n*m(2,:));
k=my*n+(mx+1);
old_to_new_dom=unique([t(4,:);k]','rows')';
old_to_new_dom=old_to_new_dom(2,:);
ee=e;
for i=1:size(old_to_new_dom,2)
    ind=(e(6,:)==i);ee(6,ind)=old_to_new_dom(i);
    ind=(e(7,:)==i);ee(7,ind)=old_to_new_dom(i);
end
tbd=(ee(7,:)==ee(6,:)); % to be deleted edges
tbk=setdiff(1:size(ee,2),find(tbd)); % to be kept edges
eee=ee(:,tbk);

tt=t;
tt(4,:)=k;