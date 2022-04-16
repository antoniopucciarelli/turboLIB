/*
TURBOMACHINERY -- LIBRARY FOR THE INITIAL TURBOMACHINERY DESIGN
AUTHOR: antonio pucciarelli 

PROGRAM DESCRIPTION
    TURBOMACHINERY 3D GEOMETRY GENERATION WITH OPENSCAD 
*/

// VARIABLES ALLOCATION
include <data.scad> 

//nRotorBlades = 1;
//nStatorBlades = 1;

// ROTOR PRINT
// variables computation
alphaRotor = 360/nRotorBlades; // computing distance angle between 2 blades
rotorStaggerAngle = 0; // setting up hub blade stagger angle
rotorOrigin = [0.001, 0, 0.97*rotorHubInletCoords[2]]; // setting up rotor hub blade section origin

// rotor generation 
for(ii = [0:nRotorBlades-1]){
    color("LightGrey", 1){
        rotate(a = ii * alphaRotor, v = [1, 0, 0]){
            translate(rotorOrigin){
                rotate(rotorStaggerAngle){
                    import(rotorName);
                };
            };
        };
    };
};

// STATOR PRINT
// variables computation
alphaStator = 360/nStatorBlades; // computing distance angle between 2 blades
statorStaggerAngle = 0; // setting up hub blade stagger angle
statorOrigin = [0.2, 0, 0.97*statorHubInletCoords[2]]; // setting up rotor hub blade section origin

// stator generation 
for(ii = [0:nStatorBlades-1]){
    color("Grey", 1.0){
        rotate(a = ii * alphaStator, v = [1, 0, 0]){
            translate(statorOrigin){
                rotate(statorStaggerAngle){
                    import(statorName);
                };
            };
        };
    };
};

// HUB PRINT
// hub generation
rotate([0,90,0]){
    color("Blue", 0.5){
        rotate_extrude(angle=360, $fn=200){
            rotate([0,0,90]){          
                polygon(points = [[0, 0], [rotorHubInletCoords[0], rotorHubInletCoords[2]], [rotorHubOutletCoords[0], rotorHubOutletCoords[2]], [statorOrigin[0], statorHubInletCoords[2]], [statorOrigin[0] + statorHubOutletCoords[0], statorHubOutletCoords[2]], [statorOrigin[0] + statorHubOutletCoords[0], 0]]);
            };
        };
    };
};