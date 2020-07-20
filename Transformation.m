%-----------------------------*ת������*----------------------------------
%-----------------------------*����ѧ*----------------------------------
function TK = Transformation(ME,Coordinates,Elength)
Coordinates_1 = Coordinates(ME(1,1),:);
Coordinates_2 = Coordinates(ME(2,1),:);
v = Coordinates_2 - Coordinates_1;
if size(Coordinates,2) == 2 %�ж��Ƕ�ά���⻹����ά����
    TK = [v(1)/Elength v(2)/Elength  0 0 ;
          0 0 v(1)/Elength v(2)/Elength];
else
    TK = [v(1)/Elength v(2)/Elength v(3)/Elength 0 0 0 ;
          0 0 0 v(1)/Elength v(2)/Elength v(3)/Elength];
end