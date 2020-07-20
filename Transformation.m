%-----------------------------*转换矩阵*----------------------------------
%-----------------------------*静力学*----------------------------------
function TK = Transformation(ME,Coordinates,Elength)
Coordinates_1 = Coordinates(ME(1,1),:);
Coordinates_2 = Coordinates(ME(2,1),:);
v = Coordinates_2 - Coordinates_1;
if size(Coordinates,2) == 2 %判断是二维问题还是三维问题
    TK = [v(1)/Elength v(2)/Elength  0 0 ;
          0 0 v(1)/Elength v(2)/Elength];
else
    TK = [v(1)/Elength v(2)/Elength v(3)/Elength 0 0 0 ;
          0 0 0 v(1)/Elength v(2)/Elength v(3)/Elength];
end