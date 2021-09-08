# Data-gap-fill
A script to fill in missing time series points in a way that preserves the characteristics of the noise in the dataset
![image](https://user-images.githubusercontent.com/39776793/132546232-6f5c83db-d535-44ee-8dcf-b70481f8d5b8.png)
# About the algorithm
The algorithm is descried in "Characterizing Frequency Stability Measurements having Multiple Data Gaps" by D. A. Howe and N. Schlossberger.
The gaps are filled via taking the nearest continuous data run of sufficient length, flipping it both horizontally and vertically, and matching the endpoints such that the resulting data is continuous.

Additional complexity is introduced in the endpoint matching.  Single-point endpoint matching can introduce problems in the presence of white noise or outliers, as demonstrated below.
![image](https://user-images.githubusercontent.com/39776793/132548378-51881758-41ad-4a5e-91a1-385c2f6048e0.png)
The result is an artificial "sawtooth" pattern in the data. To avoid this, a low pass filter is applied to each chunk of data and the endpoint matching is done according to the filtered data. The low pass takes the form
<img src="https://latex.codecogs.com/gif.latex?e%5E%7B-2%5E3%20%5Cfrac%7B%7Cs%7C%7D%7BT%7D%5Ctau_0%5E2%7D" /> 
![image](https://user-images.githubusercontent.com/39776793/132548346-c0c84436-1a38-41b5-8428-cb26b1ef1341.png)

![image](https://user-images.githubusercontent.com/39776793/132548398-38af948d-7b77-4fd7-8981-10499a49caf3.png)

# Authors
Written by Noah Schlossberger and Chloe Champagne. Chloe wrote the code to process the data, find the gaps, and address them, while Noah wrote the code to match endpoints, including the filtering function to avoid matching to outliers and high frequency noise.
