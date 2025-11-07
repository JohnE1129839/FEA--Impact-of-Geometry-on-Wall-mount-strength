%Stiffness matrix for material of uniform thickness
function ke = stiffnessMatrixUniformBar(E,h,v,n1,n2,n3)
    x1 = n1(1,1); y1 = n1(2,1);
    x2 = n2(1,1); y2 = n2(2,1);
    x3 = n3(1,1); y3 = n3(2,1);
    A = 1/2 * det([1 x1 y1; 1 x2 y2; 1 x3 y3]);
    b1 = (y2-y3)/2/A; c1 = (x3-x2)/2/A;
    b2 = (y3-y1)/2/A; c2 = (x1-x3)/2/A;
    b3 = (y1-y2)/2/A; c3 = (x2-x1)/2/A;
    B = [b1 0 b2 0 b3 0;
         0 c1 0 c2 0 c3;
         c1 b1 c2 b2 c3 b3];
    Ch = stressStrainMatrix(E,v);
    ke = h * A * (B' * Ch * B);
end



