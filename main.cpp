#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define P 12  
#define Q P
#define PI 3.142857142857
#define BATCH 320  // Number of samples per frame
#define FRAMES 5  // Number of frames for selection

int count=0; // To count the number of elements in the file

// Global arrays
double frame[FRAMES][BATCH]; // Array for storing the selected steady state frame
double Ri[FRAMES][P + 1]; // Array for storing auto correlation matrix
double Ai[FRAMES][P + 1]; // Array for storing Levinson Durbin Coefficient
double Ci[FRAMES][P + 1]; // Array for storing Cepstral Coefficients
double Cri[FRAMES][P + 1]; // Array for storing Raised Cosine Window values
double Xi[FRAMES][P + 1]; // Testing file's steady state frames
double data_arr[20000]={0.0}; // Array for storing entire data file
double Wi[P] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0}; // Tokhura Weights
double dist[FRAMES]; // Array for storing Tokhura Distances to each vowel


// For calculating DCShift
void DCshift(){
	// Calculating the mean (DC offset)
	double sum = 0;
	for (int i = 0; i < count; i++) {
		sum += data_arr[i];
	}
	double mean = sum / count;

	// Correcting DC shift by subtracting the mean from each data point
	for (int i = 0; i < count; i++) {
		data_arr[i] -= mean;
	}
}

// For doing normalization
void Normalization(){
	// Find min and max values in the data array
    double min_val = data_arr[0];
    double max_val = data_arr[0];

    for (int i = 1; i < count; i++) {
        if (data_arr[i] < min_val) {
            min_val = data_arr[i];
        }
        if (data_arr[i] > max_val) {
            max_val = data_arr[i];
        }
    }

	//Normalizing the values to the range of -5000 to 5000
	for (int i = 0; i < count; i++) {
        data_arr[i] = -5000 + (((data_arr[i] - min_val) / (max_val-min_val)) * 10000);
   }

}

// For selecting steady state frames
void SteadyStateSelection(){

	// Finding the max value's index
	double max_val = data_arr[0];
    int max_index = 0;

    for (int i = 1; i < count; i++) {
        if (data_arr[i] > max_val) {
            max_val = data_arr[i];
            max_index = i;
        }
    }

	// Taking the batch just before and just after the max value
	for (int i=0;i<5;i++){
		for (int j=0;j<320;j++){
			frame[i][j]=data_arr[max_index+j+320*(i-2)];
		}
	}
}

// Function to apply a Hamming window to the input frame
void Hamming(int i) {

    for (int n = 0; n < BATCH; n++) {

        frame[i][n] *= 0.54 - 0.46 * cos(2 * PI * n / (BATCH - 1));
    }
}

// Function to compute autocorrelation
void AutoCorrelation(int i) {

    for (int k = 0; k <= P; k++) {

		// Initialization of the first value
        Ri[i][k] = 0.0;

		// Actual calculation of Values
        for (int n = 0; n < BATCH - k; n++) {

            Ri[i][k] += frame[i][n] * frame[i][n + k];

        }
    }
}

// Function to compute LPC coefficients using Levinson-Durbin recursion
void LevinsonDurbin(int frameIdx) {
    double E[P + 1] = {0.0};        // Prediction error 
    double K[P + 1] = {0.0};        // Reflection coefficients
    double alpha[P + 1][P + 1] = {0.0}; // LPC coefficients matrix

    E[0] = Ri[frameIdx][0];         // Initial error value (E[0] is the autocorrelation at lag 0)

    for (int i = 1; i <= P; i++) {
        double sum = 0.0;

        // Compute the i-th reflection coefficient K[i]
        for (int j = 1; j < i; j++) {
            sum += alpha[j][i - 1] * Ri[frameIdx][i - j];
        }
        K[i] = (Ri[frameIdx][i] - sum) / E[i - 1];

        // Update the LPC coefficients matrix alpha
        alpha[i][i] = K[i];
        for (int j = 1; j < i; j++) {
            alpha[j][i] = alpha[j][i - 1] - K[i] * alpha[i - j][i - 1];
        }

        // Update the prediction error for the i-th order
        E[i] = (1 - K[i] * K[i]) * E[i - 1];
    }

    // Copy the LPC coefficients for the current frame to the Ai array
    for (int i = 1; i <= P; i++) {
        Ai[frameIdx][i] = alpha[i][P];
    }
}


// Function to compute Cepstral Coefficients from LPC coefficients
void ComputeCepstralCoefficients(int frameIdx) {
    // Initialize the zeroth cepstral coefficient C0
    Ci[frameIdx][0] = 0.0;
    Cri[frameIdx][0] = 0.0;

    // Compute cepstral coefficients C1 to CP
    for (int n = 1; n <= P; n++) {
        // Set the initial value of Cn to the LPC coefficient An
        Ci[frameIdx][n] = Ai[frameIdx][n];
        double sum = 0.0;

        // Calculate the contribution of previous cepstral coefficients
        for (int k = 1; k < n; k++) {
            if (n - k >= 0) {
                sum += k * Ai[frameIdx][n - k] * Ci[frameIdx][k];
            }
        }

        // Update the cepstral coefficient Cn with the accumulated sum
        Ci[frameIdx][n] += sum / n;

        // Compute the raised cosine cepstral coefficient Cri[n] for mean calculation
        Cri[frameIdx][n] += Ci[frameIdx][n] * (1 + (Q / 2) * sin(PI * n / Q));
    }
}


// Function to read signal data from file
void ReadDataFromFile(char* filename) {
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Read data from file until end-of-file or until array limit is reached
    while (!feof(file)) {
        // Read a double value from the file and store it in data_arr
        // Increment count if the value is successfully read and within bounds
        if (fscanf(file, "%lf", &data_arr[count]) == 1 && count < 20000) {
            count++;
        } else {
            break; // Exit loop if reading fails or array limit is reached
        }
    }

    // Close the file
    fclose(file);
}


// Function to write data to a file
void WriteBack(char* filename) {
    // Open the file in write mode
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Write data to the file
    for (int i = 0; i < FRAMES; i++) {
        for (int j = 1; j <= P; j++) {
            // Scale down the cepstral coefficient and write it to the file
            Cri[i][j] /= 20;
            fprintf(file, "%lf ", Cri[i][j]);
            // Reset the value in the array to 0 after writing
            Cri[i][j] = 0;
        }
        // Write a new line after each row
        fprintf(file, "\n");
    }

    // Close the file
    fclose(file);
}


// Function to read a matrix from a file into the Xi array
void RefFiles(char *filename) {
    // Open the file in read mode
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    // Read the matrix data from the file
    for (int i = 0; i < FRAMES; i++) {
        for (int j = 1; j <= P; j++) {
            // Read a double value from the file into Xi[i][j]
            if (fscanf(file, "%lf", &Xi[i][j]) != 1) {
                // If reading fails, report the error and exit
                fprintf(stderr, "Error reading data at [%d][%d]\n", i, j);
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Close the file after reading all data
    fclose(file);
}


// Function to calculate Tokhura distance between two rows
double CalcTokhura(double row1[P + 1], double row2[P + 1]) {
    double distance = 0.0;
    
    // Compute the Tokhura distance by summing weighted squared differences
    for (int i = 1; i <= P; i++) {
        double diff = row1[i] - row2[i];
        distance += Wi[i - 1] * diff * diff; // Accumulate the weighted squared difference
    }
    
    return distance;
}

// Function to compute the average Tokhura distance for all frames
double AvgTokhura() {
    double total_distance = 0.0;
    
    // Sum Tokhura distances for all frames
    for (int i = 0; i < FRAMES; i++) {
        double distance = CalcTokhura(Xi[i], Cri[i]);
        total_distance += distance; // Accumulate the total distance
    }
    
    // Compute and return the average Tokhura distance
    return total_distance / FRAMES;
}


// Function to find the index of the minimum value in the dist array
int MinIndex() {
    int x = 0; // Assume the first index (0) is the minimum initially
    
    // Loop through the dist array to find the index of the minimum value
    for (int i = 1; i < FRAMES; i++) {
        if (dist[i] < dist[x]) {
            x = i; // Update index x if a smaller value is found
        }
    }
    
    return x; // Return the index of the minimum value
}

// Function to reset the Cri array to zero
void ResetCri() {
    // Loop through each frame
    for (int i = 0; i < FRAMES; i++) {
        // Loop through each coefficient for the current frame
        for (int j = 1; j <= P; j++) {
            Cri[i][j] = 0; // Reset each coefficient to zero
        }
    }
}


int main() {

	// Define file name and output file name arrays
	char filename[100] = "";
	char outputfilename[100] = "";

	// Define an array of vowel characters for processing
	char c[5] = {'a', 'e', 'i', 'o', 'u'};

	// Loop through each vowel character
	for (int x = 0; x < 5; x++) {
		// Loop through each file number (1 to 20)
		for (int y = 1; y < 21; y++) {

			// Reinitialize global variables for each iteration
			count = 0;

			// Construct the filename for input data based on vowel and file number
			sprintf(filename, "Inputs/244101039_%c_%d.txt", c[x], y);

			// Read the signal data from the constructed filename
			ReadDataFromFile(filename);

			// Perform DC shift to remove any DC component from the signal
			DCshift();

			// Normalize the signal data to a standard range
			Normalization();

			// Select stable (steady) frames from the normalized data
			SteadyStateSelection();

			// Apply Hamming window to each frame for spectral analysis
			for (int i = 0; i < FRAMES; i++) {
				Hamming(i);
			}

			// Calculate autocorrelation, LPC coefficients, and cepstral coefficients for each frame
			for (int i = 0; i < FRAMES; i++) {
				AutoCorrelation(i);          // Compute autocorrelation for the frame
				LevinsonDurbin(i);           // Compute LPC coefficients using Levinson-Durbin
				ComputeCepstralCoefficients(i); // Compute cepstral coefficients from LPC coefficients
			}
		}

		// Construct the output filename based on the vowel character
		sprintf(outputfilename, "Outputs/244101039_%c.txt", c[x]);

		// Write the processed data to the output file
		WriteBack(outputfilename);
	}



	// Testing loop for processing files and recognizing vowels
	for (int x = 0; x < 5; x++) {
		for (int y = 1; y < 11; y++) {
			// Construct the filename for the test data based on vowel and file number
			sprintf(filename, "Tests/244101039_%c_%d.txt", c[x], y);
        
			// Reinitialize the count variable
			count = 0;

			// Read the signal data from the constructed filename
			ReadDataFromFile(filename);

			// Perform DC shift to remove any DC component from the signal
			DCshift();

			// Normalize the signal data to a standard range
			Normalization();

			// Select stable (steady) frames from the normalized data
			SteadyStateSelection();

			// Apply Hamming window to each frame for spectral analysis
			for (int i = 0; i < FRAMES; i++) {
				Hamming(i);
			}

			// Compute autocorrelation, LPC coefficients, and cepstral coefficients for each frame
			for (int i = 0; i < FRAMES; i++) {
				AutoCorrelation(i);          // Compute autocorrelation for the frame
				LevinsonDurbin(i);           // Compute LPC coefficients using Levinson-Durbin
				ComputeCepstralCoefficients(i); // Compute cepstral coefficients from LPC coefficients
			}

			// Compute Tokhura distance for each reference vowel file
			for (int i = 0; i < 5; i++) {
				// Construct the filename for the reference data
				sprintf(outputfilename, "Outputs/244101039_%c.txt", c[i]);
            
				// Read the reference data from the constructed filename
				RefFiles(outputfilename);

				// Calculate the average Tokhura distance between the test data and the reference data
				dist[i] = AvgTokhura();
			}

			// Find the index of the smallest Tokhura distance
			int z = MinIndex();
        
			// Print the recognized vowel for the current test file
			printf("%s file vowel is %c.\n", filename, c[z]);
        
			// Reset the Cri array for the next test file
			ResetCri();
		}

		// Print a new line after processing all files for the current vowel
		printf("\n");
	}

	return 0;
}

