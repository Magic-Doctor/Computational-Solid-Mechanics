%-----------------------------*转换矩阵*----------------------------------
%-----------------------------*稳定性*----------------------------------
function stabTK = stabTransformation(ME,Coordinates,Elength)
Coordinates_1 = Coordinates(ME(1,1),:);
Coordinates_2 = Coordinates(ME(2,1),:);
v = Coordinates_2 - Coordinates_1;
stabTK = [v(1)/Elength v(2)/Elength 0 0;
          -power(1-power(v(1)/Elength,2),0.5) v(1)/Elength 0 0;
          0 0 v(1)/Elength v(2)/Elength;
          0 0 -power(1-power(v(1)/Elength,2),0.5) v(1)/Elength];
%平面稳定性问题的转换矩阵(几何/初应力矩阵)为
%[cos(x,x') cos(y,x') 0 0;
% -sin(x,x')  cos(x,x') 0 0;
% 0 0 cos(x,x') cos(y,x') ;
% 0 0 -sin(x,x')  cos(x,x')]
% 其中sin(x,x') = -根号下(1-cos(x,x')平方);