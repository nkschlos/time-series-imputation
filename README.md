# Data-gap-fill
A script to fill in missing time series points in a way that preserves the characteristics of the noise in the dataset
![image](https://user-images.githubusercontent.com/39776793/132546232-6f5c83db-d535-44ee-8dcf-b70481f8d5b8.png)
# How to use the script

Ensure that fillgaps.py, runFillgaps.py, and the csv or txt file you want to fill are in the same folder. 
Also make sure that the csv or text file is formatted as two (comma delimited) columns, the first being the time of the data collected in MJD, and the second being the data, nominally the clock residual in microseconds.

Simply change the name of the file in runFillgaps.py and run this file in python.

# About the algorithm
The algorithm is descried in "Characterizing Frequency Stability Measurements having Multiple Data Gaps" by D. A. Howe and N. Schlossberger.
The gaps are filled via taking the nearest continuous data run of sufficient length, flipping it both horizontally and vertically, and matching the endpoints such that the resulting data is continuous.

Additional complexity is introduced in the endpoint matching.  Single-point endpoint matching can introduce problems in the presence of white noise or outliers, as demonstrated below.

![image](https://user-images.githubusercontent.com/39776793/132553326-19be783c-da27-4445-8b4c-98b9e859d055.png | width=100)

The result is an artificial "sawtooth" pattern in the data. To avoid this, a low pass filter is applied to each chunk of data and the endpoint matching is done according to the filtered data. The low pass takes the form

e^(-(2^3)(tau_0^2)|f|/T)

where tau_0 is the time sampling period, T is the length of the data being filtered, and f is the discrete frequency variable.

![image](https://user-images.githubusercontent.com/39776793/132553213-358285a6-65dd-4c16-8ad0-40e3470b984a.png | width=100)


A final complication arises because the fft treats the data as periodic, so it will try to match the first point to the last, which can result in an up-tick or down-tick at the end of the data, where we are trying to match. The solution is to pad the data with a flipped copy of itself on either side before filtering, and use the points of reflection for the matching.

![image](https://user-images.githubusercontent.com/39776793/132553286-23cb80fa-73fa-4a8d-87dc-11b23307b7e7.png | width=100)

# Authors
Written by Noah Schlossberger and Chloe Champagne. Chloe wrote the code to process the data, find the gaps, and address them, while Noah wrote the code to match endpoints, including the filtering function to avoid matching to outliers and high frequency noise.
