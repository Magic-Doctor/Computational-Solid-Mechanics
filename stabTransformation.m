%-----------------------------*ת������*----------------------------------
%-----------------------------*�ȶ���*----------------------------------
function stabTK = stabTransformation(ME,Coordinates,Elength)
Coordinates_1 = Coordinates(ME(1,1),:);
Coordinates_2 = Coordinates(ME(2,1),:);
v = Coordinates_2 - Coordinates_1;
stabTK = [v(1)/Elength v(2)/Elength 0 0;
          -power(1-power(v(1)/Elength,2),0.5) v(1)/Elength 0 0;
          0 0 v(1)/Elength v(2)/Elength;
          0 0 -power(1-power(v(1)/Elength,2),0.5) v(1)/Elength];
%ƽ���ȶ��������ת������(����/��Ӧ������)Ϊ
%[cos(x,x') cos(y,x') 0 0;
% -sin(x,x')  cos(x,x') 0 0;
% 0 0 cos(x,x') cos(y,x') ;
% 0 0 -sin(x,x')  cos(x,x')]
% ����sin(x,x') = -������(1-cos(x,x')ƽ��);