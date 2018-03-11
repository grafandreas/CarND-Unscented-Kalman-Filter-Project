# Unscented Kalman Filter Project

This is the code for the unscented Kalman filter project of the SDCND. The structure of this file
follows the rubric.

## Compilation

No changes were made to the build process, the code should compile out-of-the box on all supported platforms.

## Accuracy


When running with the simulator, the RMSE shown is as follows:

Variable|Value
--------|----------
x|0.0736 
y|0.0829
vx|0.3423
vx|0.2590

## Code Structure

The basic code structure of the initial code has been slightly modified to be somewhat more modular. The code has dedicated 
methods for:

* Initialization : void UKF::Init(MeasurementPackage measurement_pack)  
* Prediciting Sigma Points: UKF::PredictSigmaPoints(double delta_t)
* Updating the Covariance Matrix : UKF::PredictSigmaPoints(double delta_t)

In addition, the common code for updating Lidar And radar has been factored out to UKF::UpdateUKF

Note that for deployment in a real embedded device, some of these refactorings could be unrolled, to avoid the overhead of function calls. This would probably been done by using the "inline" keyword to still keep the code modular.

### Optimizations
Some code has been refactored, especially calculation of power of 2, so that the value has to be calculated only once. Tigonometric calculations have also been extracted, assuming that the performance of these calls can be quite slow.

### Debug code
Assertions have been added to detect NaNs early. These are usually not problematic, since they can be deactivated by compiler switchesfor production code.

 