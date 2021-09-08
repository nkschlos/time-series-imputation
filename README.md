# Data-gap-fill
A script to fill in missing time series points in a way that preserves the characteristics of the noise in the dataset
![image](https://user-images.githubusercontent.com/39776793/132546232-6f5c83db-d535-44ee-8dcf-b70481f8d5b8.png)
# About the algorithm
The algorithm is descried in "Characterizing Frequency Stability Measurements having Multiple Data Gaps" by D. A. Howe and N. Schlossberger.
The gaps are filled via taking the nearest continuous data run of sufficient length, flipping it both horizontally and vertically, and matching the endpoints such that the resulting data is continuous.
# Authors
Written by Noah Schlossberger and Chloe Champagne. Chloe wrote the code to process the data, find the gaps, and address them, while Noah wrote the code to match endpoints, including the filtering function to avoid matching to outliers and high frequency noise.
