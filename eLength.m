%-----------------------------*单元长度*----------------------------------
function Elength = eLength(ME,Coordinates)
NE = length(ME);
Elength = zeros(1,NE);
for i = 1:NE
    Coordinates1 = Coordinates(ME(1,i),:);
    Coordinates2 = Coordinates(ME(2,i),:);
    Elength(i) = sqrt(sum(power(Coordinates1 - Coordinates2,2)));
end