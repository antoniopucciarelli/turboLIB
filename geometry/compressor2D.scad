include <rotorCoords.scad>
include <statorCoords.scad>

rotorOrigin1       = [0, 0, 0];
nRotorBlades1      = 5; 
rotorPitch1        = 0.75;
rotorStaggerAngle1 = 30;

for(ii = [0:nRotorBlades1]){
    color("LightGrey", 1.0){
        translate(rotorOrigin1 + [0, ii * rotorPitch1, 0]){
        linear_extrude(height = 2)
            rotate(rotorStaggerAngle1)
            polygon(points = rotor1);
        };
    };
};

statorOrigin1         = [3, 0, -0.5];
nStatorBlades1        = 5;
statorPitch1          = 1; 
statorStaggerAngle1   = -60;
rotorStatorClearance1 = 1;

for(ii = [0:nStatorBlades1]){
    color("SlateGrey", 1.0){
        translate(statorOrigin1 + [0, ii * statorPitch1, 0]){
        linear_extrude(height = 3)
            rotate(statorStaggerAngle1)
            polygon(points = stator1);
        };
    };
};