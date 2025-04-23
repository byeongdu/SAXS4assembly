function calcInt = calculate_3D_polynomial_Rogers(X,Y,Z,Pos,Order,Coeff)
%https://github.com/AET-MetallicGlass/Supplementary-Data-Codes/blob/master/4_Tracing_and_classification/src/calculate_3D_polynomial_Rogers.m
    calcInt = zeros(size(X));
    for i=1:size(Order,1)
       calcInt = calcInt + Coeff(i)*(X-Pos(1)).^Order(i,1).*(Y-Pos(2)).^Order(i,2).*(Z-Pos(3)).^Order(i,3);
    end
end