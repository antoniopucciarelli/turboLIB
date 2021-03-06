# turboLIB

**Python library** for the **preliminary design** of a turbomachinery stage.

The library is subdivided into different modules that allow to design the stage. 

## Modules
- ``` turboClass ``` : this class generates the blade and the main thermodynamics quantities along the blade span. 
    
    * ``` blade ``` is the blade object and stores all the operations needed for the study of the **NISRE**.

- ``` turboCoeff ``` this class stores all the modules needed for:
    
    - **shape optimization**: ``` leiblein ```
    - **losses computation**: ``` losses ```
    - **adimensional desing**: ``` similarity ```  
    
- ``` geometry ``` this class generates the actual blade shape
    - **airfoil generator**: ```geometryData```

## Compressor design

The initial design steps are related to the choice of the **number of blades**, this choice has been made using ```bladeDesign.py``` outputs. 

The **mean line** study is made with ```similarityStudy.py``` script and it allows to set the **mean flow behaviour** to follow for the generation of the **3D blades**.

The file ```compressorDesign.py``` uses the ```turboLIB``` program for the preliminary design of a compressor.

In the ```compressorDesign.py``` file there are the **initial conditions** and **constraints** of the compressor. 

```compressorDesign.py``` will save the output text in ```compressor_<rD>_<rMean>_<nRotorBlades>_<nStatorBlades>.txt``` and the blade geometry in ```.stl``` format into ```container/```. 

At the end a ```.scad``` file is generated and it can be used with ```openSCAD```.

A summary of the preliminary compressor design is explained by [these slides](https://github.com/antoniopucciarelli/turboLIB/blob/main/latex/main.pdf).
