# Vowel Recognition Using Cepstral Coefficients and Tokhura's Distance

## Overview
This project processes audio files of vowels, extracts cepstral coefficients, and predicts the vowel spoken using Tokhura’s distance. The system performs key signal processing steps like DC Shift correction, normalization, steady-state selection, Hamming windowing, autocorrelation, LPC coefficient calculation, cepstral coefficient extraction, and raised sine window application.

## Features
1. **Input Files**: The project processes text files representing audio signals of vowels.
2. **DC Shift & Normalization**: Corrects for DC bias and normalizes the input audio.
3. **Steady State Selection**: Extracts steady-state portions of the vowel sounds.
4. **Hamming Window**: Reduces spectral leakage in the signal.
5. **Autocorrelation & LPC Coefficients**: Computes LPC coefficients using autocorrelation.
6. **Cepstral Coefficients**: Extracts cepstral coefficients for further analysis.
7. **Vowel Prediction**: Utilizes **Tokhura's distance** to predict the vowel from test files.

## Project Structure
- `main.cpp`: The core C++ file that implements the vowel recognition system.
- `input/`: Folder containing input text files of vowel sounds.
- `output/`: Folder where the processed results (e.g., cepstral coefficients) are saved.
- `test/`: Folder with test files used for vowel prediction.

## Prerequisites
- **C++ Compiler**: Ensure you have a C++ compiler installed (e.g., GCC, Visual Studio).

## How to Run
1. **Clone or Download** the project files to your local machine.
2. **Compile** the `main.cpp` file using any C++ compiler:
   ```bash
   g++ main.cpp -o vowel_recognition
   ```
3. **Run the compiled program**:
   ```bash
   ./vowel_recognition
   ```
   - The program will process files in the `input/` folder, generate outputs in the `output/` folder, and predict vowels from files in the `test/` folder.
4. **View Results**: The predicted vowel and other results will be displayed in the terminal, and detailed outputs will be saved in the `output/` folder.

## Key Algorithms
- **DC Shift and Normalization**: Adjusts the audio signals for consistency.
- **Steady State Selection**: Extracts the most stable segment of the vowel sound.
- **Hamming Window**: Applied to minimize spectral leakage.
- **Autocorrelation & LPC Coefficients**: Calculates Linear Predictive Coding coefficients.
- **Cepstral Coefficients**: Extracted from LPC coefficients for effective vowel classification.
- **Tokhura's Distance**: Measures the similarity between test and reference vowels to make predictions.

## Results
- The results will be displayed in the terminal, showing the predicted vowel for each test file.
- Processed data (like cepstral coefficients) will be saved in the `output/` folder.

## License
This project is licensed under the MIT License.
