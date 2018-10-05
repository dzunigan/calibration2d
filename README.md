# Automatic Multisensor Extrinsic Calibration from Egomotion
Implementation of an automatic method to estimate the 2D extrinsic calibration parameters of multiple sensors through their observed egomotion, explicitly handling scale ambiguity from cameras.

**Authors:** [David Zuñiga-Noël](http://mapir.isa.uma.es/mapirwebsite/index.php/people/270), [Jose-Raul Ruiz-Sarmiento](http://mapir.uma.es/mapirwebsite/index.php/people/108-jose-raul-ruiz-sarmiento), [Ruben Gomez-Ojeda](http://mapir.isa.uma.es/mapirwebsite/index.php/people/164-ruben-gomez), and [Javier Gonzalez-Jimenez](http://mapir.isa.uma.es/mapirwebsite/index.php/people/95-javier-gonzalez-jimenez)

**License:**  [GPLv3](https://raw.githubusercontent.com/dzunigan/calibration2d/master/LICENSE.txt)

## 1. Dependencies

* Boost (1.58.0.1ubuntu1)
   ```
   sudo apt install libboost-all-dev
   ```
* Ceres (1.14.0-facb199)

   See ceres [documentation](http://ceres-solver.org/installation.html#linux).
   
* CMake (3.5.1-1ubuntu1)
   ```
   sudo apt install cmake
   ```
* Eigen (3.3~beta1-2)
   ```
   sudo apt install libeigen3-dev
   ```
* Gflags (2.1.2-3)
   ```
   sudo apt install libgflags-dev
   ```
* Glog (0.3.4-0.1)
   ```
   sudo apt install libgoogle-glog-dev
   ```

## 2. Data Format

The input are synchronized per-sensor incremental 2D motions. For each sensor, incremental motions are grouped into a single text file, where each line contains a synchronized 2D incremental motion (angles in radians) in the following format:
```
x y yaw
```

Note: all sensors must have the same number of observations

To compute synchronized incremental motions from estimated trajectories, please see [tools/sync](tools/sync).

## 3. Usage

The multisensor calibration method can be invoked as:
```
calibrate [options] <motions1> <motions2> ...
```

Notes:

* The method calibrates all sensors with respect to a reference one. The reference sensor is assumed to be the one observing the first incremetal motions.

* Use the `--scale_ambiguous=n,m,...` option to indicate which motions have scale ambiguity.

## 4. Example

**TODO**
