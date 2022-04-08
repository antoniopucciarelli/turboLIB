statorOrigin1 = [0.1, 0, 0.2];

// HUB
// variable allocation
rMean              = 0.325;
// rotor properties
rotorHubChord      = 0.0627869302309185;
rotorHub_xInlet    = 0.0;
rotorHub_yInlet    = 0.0;
rotorHub_zInlet    = 0.20086511551156916;
rotorHub_xOutlet   = 0.0627869302309185;
rotorHub_yOutlet   = 0.014470943474431435;
rotorHub_zOutlet   = 0.2483654691199017;
rotorHeigyhInlet   = 0.2482697689768617;
rotorHeigthOutlet  = 0.19907458767017935;
// stator properties
statorHubChord      = 0.07446010544196538;
statorHub_xInlet    = 0.0;
statorHub_yInlet    = 0.0;
statorHub_zInlet    = 0.22546270616491032;
statorHub_xOutlet   = 0.07446010544196538;
statorHub_yOutlet   = 0.030290561214741725;
statorHub_zOutlet   = 0.24235713251146196;
statorHeigthInlet   = 0.19907458767017935;
statorHeigthOutlet  = 0.18212580449864338;

// hub shape generation 
rotate([0,90,0]){
    color("Blue", 0.5){
        rotate_extrude(angle=360, $fn=200){
            rotate([0,0,90]){          
                polygon(points = [[0, 0], [rotorHub_xInlet, rotorHub_zInlet], [rotorHub_xOutlet, rotorHub_zOutlet], [statorOrigin1[0], statorHub_zInlet], [statorOrigin1[0] + statorHub_xOutlet, statorHub_zOutlet], [statorOrigin1[0] + statorHub_xOutlet, 0]]);
            };
        };
    };
};