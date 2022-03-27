include <rotorCoords.scad>
include <statorCoords.scad>

nRotorBlades1 = 30; 
alphaRotor1 = floor(360/nRotorBlades1);
rotorStaggerAngle1 = 30;
rotorOrigin1 = [0, 0, 4];
rotorBladeHeight1 = 2;

for(ii = [0:nRotorBlades1-1]){
    color("LightGrey", 1.0){
        rotate(a = ii * alphaRotor1, v = [1, 0, 0]){
            translate(rotorOrigin1){
                linear_extrude(height=rotorBladeHeight1)
                    rotate(rotorStaggerAngle1)
                    polygon(points = rotor1);
            };
        };
    };
};

nStatorBlades1 = 30; 
alphaStator1 = floor(360/nStatorBlades1);
statorStaggerAngle1 = -60;
statorOrigin1 = [2, 0, 4];
statorBladeHeight1 = 2;

for(ii = [0:nStatorBlades1-1]){
    color("SlateGrey", 1.0){
        rotate(a = ii * alphaStator1, v = [1, 0, 0]){
            translate(statorOrigin1){
                linear_extrude(height=statorBladeHeight1)
                    rotate(statorStaggerAngle1)
                    polygon(points = stator1);
            };
        };
    };
};