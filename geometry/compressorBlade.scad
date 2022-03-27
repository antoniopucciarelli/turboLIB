nRotorBlades1 = 20; 
alphaRotor1 = floor(360/nRotorBlades1);
rotorStaggerAngle1 = 10;
rotorOrigin1 = [0, 0, 4];
rotorBladeHeight1 = 2;

for(ii = [0:nRotorBlades1-1]){
    color("LightGrey", 1.0){
        rotate(a = ii * alphaRotor1, v = [1, 0, 0]){
            translate(rotorOrigin1){
                rotate(rotorStaggerAngle1){
                    import("cad.stl");
                };
            };
        };
    };
};
