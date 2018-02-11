# ConfFactors
![alt text](https://github.com/blupiac/ConfFactors/raw/master/report/mosaicdillo.png)  
  
Implementation of the paper: [Characterizing Shape Using Conformal Factors](http://www.cs.technion.ac.il/~gotsman/AmendedPubl/Miri/Shape08.pdf)  
  
Read the report of the project [here](https://github.com/blupiac/ConfFactors/blob/master/report/IG3DA_Report.pdf).  

## Using the application

### Inputs

The application can be built and executed without any arguments, in this case, the default mesh (*ie.* armadillo1) will be used for the calculations. The path to a .off can also be given as an argument, this file is then used for the calculations.

### Outputs

1. A visual output, in the shape of an OpenGL window, that displays the mesh with the conformal factor as vertex colors. The same color scale as in the original paper is used. It can be replaced by a greyscale by commenting the pre-processor constant "GRADIENT" in Main.cpp, l.32. The mesh can be rotated by moving the mouse with the left button pressed, panned by moving the mouse with the right button pressed and zoomed in/out using the mouse wheel.  
    
2. A .csv file containing the signature of the mesh in the output folder: this file can be used as as input for different applications that use the conformal factor.


### Additional features

- By turning on the TIMER macro (mesh.cpp l.22, mesh.h l.29), the computation time of different parts of the application are logged to the console.  
- By turning on the DEBUG macro (mesh.cpp l.21, mesh.h l.28, Main.cpp l.33), the visual output is the gaussian curvature instead of the conformal factor.  
- The PRUNEFACT macro (mesh.cpp l.28) determines the percentage of values pruned from the conformal factors vector. Using 0 disabled pruning, and 0.01 for example removes 0.5% of the highest and lowest values from the vector.  
- The NOISEFACT macro (mesh.cpp l.29) determines the intensity of noise added to the mesh. Try values between 0.01 and 1.  
