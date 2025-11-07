%Stress-strain Matrix used for a uniform material
function C = stressStrainMatrix(E,v)
    C = [E/(1-v^2), v*E/(1-v^2), 0; v*E/(1-v^2), E/(1-v^2), 0; 0,0, E/(2+2*v)];
end