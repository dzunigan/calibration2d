# Automatic Multisensor Extrinsic Calibration from Egomotion
Implementation of an automatic method to estimate the 2D extrinsic calibration parameters of multiple sensors through their observed egomotion, explicitly handling scale ambiguity from cameras.

**Authors:** [David Zuñiga-Noël](http://mapir.isa.uma.es/mapirwebsite/index.php/people/270), [Jose-Raul Ruiz-Sarmiento](http://mapir.uma.es/mapirwebsite/index.php/people/108-jose-raul-ruiz-sarmiento), [Ruben Gomez-Ojeda](http://mapir.isa.uma.es/mapirwebsite/index.php/people/164-ruben-gomez), and [Javier Gonzalez-Jimenez](http://mapir.isa.uma.es/mapirwebsite/index.php/people/95-javier-gonzalez-jimenez)

**License:**  [GPL v3](https://raw.githubusercontent.com/dzunigan/calibration2d/master/LICENSE.txt)

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

## 3. Usage

