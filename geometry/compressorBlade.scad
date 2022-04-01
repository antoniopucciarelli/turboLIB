/*
TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
AUTHOR: antonio pucciarelli 

PROGRAM DESCRIPTION
    TURBOMACHINERY 3D GEOMETRY GENERATION WITH OPENSCAD 
*/

// variables allocation 
nRotorBlades1 = 40;                         // number of blades
alphaRotor1 = floor(360/nRotorBlades1);     // computing distance angle between 2 blades
rotorStaggerAngle1 = 10;                    // setting up hub blade stagger angle
rotorOrigin1 = [0, 0, 4];                   // setting up rotor hub blade section origin

// rotor generation 
for(ii = [0:nRotorBlades1-1]){
    color("LightGrey", 1.0){
        rotate(a = ii * alphaRotor1, v = [1, 0, 0]){
            translate(rotorOrigin1){
                rotate(rotorStaggerAngle1){
                    import("../container/cad.stl");
                };
            };
        };
    };
};
