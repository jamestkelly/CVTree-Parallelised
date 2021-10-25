/// ------------------------------------------------------------------------------------------------------------- ///
///											CVTree [Improved] | Parallelised									  ///
/// ------------------------------------------------------------------------------------------------------------- ///
/// Author: Jim Tran Kelly | N9763686																			  ///
/// ------------------------------------------------------------------------------------------------------------- ///
/// Description: This script has been updated and parallelised for the purposes of CAB401 High Performance and    ///
/// Parallel Computing at the Queensland University of Technology (QUT), Semester 2, 2021. The script in its      ///
/// original form is a program to take a file as input through the command-line, containing the number of strings ///
/// in the first line and following that several bacteria names. These names are then concatenated with .faa	  ///
/// matching to files in the /data/ directory with gene names and their respective gene represented as DNA. The   ///
/// program collects the "kmer" or subset of an entire gene and then compares that against other bacteria "kmers" ///
/// and finds the correlation between any given bacteria. This application was found on the Blackboard Assignment ///
/// tab for CAB401 as one of the provided projects that is available for parallelization. This program has been   ///
/// parallelised using OpenMP so as to be able to split for-loop operations across multiple threads or cores.     ///
/// Additionally, there have been several changes made to the provided source code so as to improve the best      ///
/// sequential version of the application.																		  ///
/// ------------------------------------------------------------------------------------------------------------- ///
/// Version: 10.2																								  ///
/// ------------------------------------------------------------------------------------------------------------- ///
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory.h>
#include <chrono>
#include <omp.h>
using namespace std;

// Declare global variables
int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

/// <summary>
/// Converts a two dimensional (2D) vector to a one dimensional (1D) form for simplified data processing. This
/// method iterates through all elements in the supplied 2D vector and appends them to a 1D vector for output.
/// </summary>
/// <param name="corrVector">
/// A two dimensional (2D) vector of doubles.
/// </param>
/// <returns>
/// A one dimensional (1D) vector of doubles.
/// </returns>
vector<double> MakeVectorOneDimension(vector<vector<double>> corrVector) {
	vector<double> resultVector; // Initialise one dimensional vector of doubles

	// Iterate through all values in the input vector
	for (int i = 0; i < corrVector.size(); i++) {
		for (int j = 0; j < corrVector[i].size(); j++) {
			if (corrVector[i][j] != 0.0) { // Ignore 0.0 values
				resultVector.push_back(corrVector[i][j]); // Append to the output vector
			}
		}
	}

	return resultVector; // Return the resulting one dimensional vector of doubles
}

/// <summary>
/// Method to write a two dimensional (2D) vector of doubles to the specified file. This method overwrites
/// the original contents of the file, or creates a new one with the same name, depending on the presence of
/// a file matching in name to the parameter "fileName".
/// </summary>
/// <param name="corrArr">
/// The two dimensional (2D) vector of doubles to be written to a file on disk.
/// </param>
/// <param name="fileName">
/// The file name or path to be written to.
/// </param>
void WriteToFile(vector<vector<double>> corrVector, string fileName) {
	std::ofstream outputFile(fileName); // Set output file

	// Convert the 2D vector to 1D
	vector<double> outputVector = MakeVectorOneDimension(corrVector);

	// Iterate through all values and write them to the file
	for (const auto& value : outputVector) {
		outputFile << std::defaultfloat
			<< std::setprecision(std::numeric_limits<long double>::digits10)
			<< value << "\n";
	}
}

/// <summary>
/// Method to read in a vector of doubles from a given file and store it for verification of parallel
/// results.
/// </summary>
/// <param name="fileName">
/// The file name or path to be read from.
/// </param>
/// <returns>
/// A one dimensional (1D) vector of doubles.
/// </returns>
vector<double> ReadVectorFromFile(string fileName) {
	double num = 0.0; // Initialise variable to store line contents
	vector<double> inputResult; // Initialse a 1D vector to store file data

	// Set input file
	std::ifstream inputFile(fileName, std::ios::in);

	// Check to see the file was opened correctly
	if (!inputFile.is_open()) {
		std::cerr << "There was a problem opening the input file.\n"; // Print error message
		exit(1); // Exit
	}

	// Store values from input file until end of file (EOF)
	while (inputFile >> num) {
		inputResult.push_back(num);
	}

	return inputResult; // Return the resulting vector of doubles
}

/// <summary>
/// Method to compare two vectors of results against each other. This method is used to verify that the
/// parallelised code's results match to those of the sequential program.
/// </summary>
/// <param name="seqResult">
/// A one dimensional (1D) vector containing doubles representing the sequential program's results.
/// </param>
/// <param name="parResult">
/// A one dimensional (1D) vector containing doubles representing the parallel program's results.
/// </param>
void CompareResults(vector<double> seqResult, vector<double> parResult) {
	double difference = 0.0; // Initialise difference variable to 0
	int mismatchCount = 0; // Initialise a counter for the number of mismatches

	// Iterate through all pairs of values
	for (int i = 0; i < seqResult.size(); i++) {
		difference = abs(seqResult[i] - parResult[i]); // Calculate the difference

		// If there is a significant difference in correlation
		if (difference >= EPSILON) {
			printf("Mismatch of: %.20f at position %d\n", difference, i);
			printf("Sequential: %.20f | Parallel: %.20f\n", seqResult[i], parResult[i]);
			mismatchCount++; // Increment the mismatches encountered
		}
	}

	// Print corresponding message to the number of mismatches encountered
	if (mismatchCount == 0) {
		printf("The two vectors match, no mismatches encountered.\n");
	}
	else {
		printf("A total of %d mismatches between the sequential and parallel results were found.\n", mismatchCount);
	}
}

/// <summary>
/// Method to set the number of threads for the program to use.
/// </summary>
/// <param name="numThreads">
/// An integer indicating the number of threads to use.
/// </param>
void SetThreads(int numThreads) {
	omp_set_num_threads(numThreads); // Set the number of threads
	printf("Running with %d threads.\n", numThreads);
}

/// <summary>
/// Method to initialise variabels.
/// </summary>
void Init()
{
	M2 = 1;
	for (int i = 0; i < LEN - 2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
		M2 *= AA_NUMBER;
	M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
	M = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
}

/// <summary>
/// Class for the instantiation of 'Bacteria' objects.
/// </summary>
class Bacteria
{
private:
	long* vector;
	long* second;
	long one_l[AA_NUMBER];
	long indexs;
	long total;
	long total_l;
	long complement;

	/// <summary>
	/// Method to initialise vectors for creation of 'Bacteria' objects.
	/// </summary>
	void InitVectors()
	{
		vector = new long[M];
		second = new long[M1];
		memset(vector, 0, M * sizeof(long));
		memset(second, 0, M1 * sizeof(long));
		memset(one_l, 0, AA_NUMBER * sizeof(long));
		total = 0;
		total_l = 0;
		complement = 0;
	}

	/// <summary>
	/// Method to initialse a buffer.
	/// </summary>
	/// <param name="buffer">
	/// An array of 'chars'.
	/// </param>
	void init_buffer(char* buffer)
	{
		complement++;
		indexs = 0;
		for (int i = 0; i < LEN - 1; i++)
		{
			short enc = encode(buffer[i]);
			one_l[enc]++;
			total_l++;
			indexs = indexs * AA_NUMBER + enc;
		}
		second[indexs]++;
	}

	/// <summary>
	/// Method to continue the buffer.
	/// </summary>
	/// <param name="ch">
	/// A character or 'char'.
	/// </param>
	void cont_buffer(char ch)
	{
		short enc = encode(ch);
		one_l[enc]++;
		total_l++;
		long index = indexs * AA_NUMBER + enc;
		vector[index]++;
		total++;
		indexs = (indexs % M2) * AA_NUMBER + enc;
		second[indexs]++;
	}

public:
	long count;
	std::vector<double> tv; // Initiallise tv as a vector of doubles
	std::vector<long> ti; // Initialise ti as a vector of longs

	/// <summary>
	/// Constructor method for the creation of 'Bacteria' objects.
	/// </summary>
	/// <param name="filename">
	/// An array of 'char', containing the filename of the given 'Bacteria'.
	/// </param>
	Bacteria(char* filename)
	{
		FILE* bacteria_file;
		errno_t OK = fopen_s(&bacteria_file, filename, "r");

		if (OK != 0)
		{
			fprintf(stderr, "Error: failed to open file %s\n", filename);
			exit(1);
		}

		InitVectors();

		char ch;
		while ((ch = fgetc(bacteria_file)) != EOF)
		{
			if (ch == '>')
			{
				while (fgetc(bacteria_file) != '\n'); // skip rest of line

				char buffer[LEN - 1];
				fread(buffer, sizeof(char), LEN - 1, bacteria_file);
				init_buffer(buffer);
			}
			else if (ch != '\n')
				cont_buffer(ch);
		}

		long total_plus_complement = total + complement;
		double total_div_2 = total * 0.5;
		int i_mod_aa_number = 0;
		int i_div_aa_number = 0;
		long i_mod_M1 = 0;
		long i_div_M1 = 0;

		double one_l_div_total[AA_NUMBER];
		for (int i = 0; i < AA_NUMBER; i++)
			one_l_div_total[i] = (double)one_l[i] / total_l;

		double* second_div_total = new double[M1];
		for (int i = 0; i < M1; i++)
			second_div_total[i] = (double)second[i] / total_plus_complement;

		count = 0;

		// Reserve memory for the vectors
		tv.reserve(M);
		ti.reserve(M);

		for (long i = 0; i < M; i++)
		{
			double p1 = second_div_total[i_div_aa_number];
			double p2 = one_l_div_total[i_mod_aa_number];
			double p3 = second_div_total[i_mod_M1];
			double p4 = one_l_div_total[i_div_M1];
			double stochastic = (p1 * p2 + p3 * p4) * total_div_2;

			if (i_mod_aa_number == AA_NUMBER - 1)
			{
				i_mod_aa_number = 0;
				i_div_aa_number++;
			}
			else
				i_mod_aa_number++;

			if (i_mod_M1 == M1 - 1)
			{
				i_mod_M1 = 0;
				i_div_M1++;
			}
			else
				i_mod_M1++;

			if (stochastic > EPSILON)
			{
				tv.push_back((vector[i] - stochastic) / stochastic);
				ti.push_back(i);
				count++;
			}
		}

		// Shrink the vectors to the number of elements
		tv.shrink_to_fit();
		ti.shrink_to_fit();

		// 'Delete' the vectors & other data structures
		delete[] second_div_total; // Changed to delete[] to remove Warning C6283
		delete vector;
		delete second;
		tv.clear();
		ti.clear();

		fclose(bacteria_file);
	}
};

/// <summary>
/// Method to read in a file containing information about a 'Bacteria'.
/// </summary>
/// <param name="input_name">
/// The input file name.
///</param>
void ReadInputFile(const char* input_name)
{
	FILE* input_file;
	errno_t OK = fopen_s(&input_file, input_name, "r");

	if (OK != 0)
	{
		fprintf(stderr, "Error: failed to open file %s (Hint: check your working directory)\n", input_name);
		exit(1);
	}

	fscanf_s(input_file, "%d", &number_bacteria);
	bacteria_name = new char* [number_bacteria];

	for (long i = 0; i < number_bacteria; i++)
	{
		char name[10];
		fscanf_s(input_file, "%s", name, 10);
		bacteria_name[i] = new char[20];
		sprintf_s(bacteria_name[i], 20, "data/%s.faa", name);
	}
	fclose(input_file);
}

/// <summary>
/// Method to compare two given 'Bacteria' objects. This method calculates the correlation between
/// the two 'Bacteria' objects.
/// </summary>
/// <param name="b1">
/// A 'Bacteria' object.
/// </param>
/// <param name="b2">
/// A 'Bacteria' object.
/// </param>
/// <returns>
/// A correlation score.
/// </returns>
double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
	double correlation = 0;
	double vector_len1 = 0;
	double vector_len2 = 0;
	long p1 = 0;
	long p2 = 0;
	while (p1 < b1->count && p2 < b2->count)
	{
		long n1 = b1->ti[p1];
		long n2 = b2->ti[p2];
		if (n1 < n2)
		{
			double t1 = b1->tv[p1];
			vector_len1 += (t1 * t1);
			p1++;
		}
		else if (n2 < n1)
		{
			double t2 = b2->tv[p2];
			p2++;
			vector_len2 += (t2 * t2);
		}
		else
		{
			double t1 = b1->tv[p1++];
			double t2 = b2->tv[p2++];
			vector_len1 += (t1 * t1);
			vector_len2 += (t2 * t2);
			correlation += t1 * t2;
		}
	}
	while (p1 < b1->count)
	{
		long n1 = b1->ti[p1];
		double t1 = b1->tv[p1++];
		vector_len1 += (t1 * t1);
	}
	while (p2 < b2->count)
	{
		long n2 = b2->ti[p2];
		double t2 = b2->tv[p2++];
		vector_len2 += (t2 * t2);
	}

	return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

/// <summary>
///	Method to compare all 'Bacteria' objects as contained within the ./data/ directory. This method calls
/// CompareBacteria() and determines the correlation between any two given 'Bacteria' objects. Additionally,
/// using OpenMP the methods loads the 'Bacteria' objects and compares them using parallel processing.
/// </summary>
/// <returns>
/// A two dimensional (2D) vector of doubles containing the correlation values between the given 'Bacteria'
/// objects.
/// </returns>
vector<vector<double>> CompareAllBacteria() {
	// Initialise and allocate a vector
	vector<vector<double>> corrVector(number_bacteria);

	// Iterate through the number of bacteria and initialise a 'row' to store values
	for (int i = 0; i < number_bacteria; i++) {
		corrVector[i] = vector<double>(number_bacteria, 0.0);
	}

	// Initialise two vectors of integers to store indices
	vector<int> I;
	vector<int> J;

	// Populate the vectors
	for (int i = 0; i < number_bacteria; i++) {
		for (int j = i + 1; j < number_bacteria; j++) {
			I.push_back(i); // Append the index 'i' to the vector
			J.push_back(j); // Append the index 'j' to the vector
		}
	}

	// Instantiate array of pointers to 'Bacteria' objects
	Bacteria** b = new Bacteria * [number_bacteria];

#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		// Iterate through all bacteria and load them
		for (int i = 0; i < number_bacteria; i++)
		{
			b[i] = new Bacteria(bacteria_name[i]);
		}

#pragma omp for schedule(dynamic)
		// Using a one-dimensional iteration space
		for (int k = 0; k < I.size(); k++) {
			int i = I[k]; // Fetch the 'i' index
			int j = J[k]; // Fetch the 'j' index
			corrVector[i][j] = CompareBacteria(b[i], b[j]); // Calculate and store correlation value
		}

	}

	return corrVector; // Return the resulting 2D vector
}

/// <summary>
/// Method to compare all 'Bacteria' objects as contained within the ./data/ directory. This method calls
/// CompareBacteria() and determines the correlation between any two given 'Bacteria' objects. Notably,
/// this method performs this comparison sequentially.
/// </summary>
/// <returns>
/// A two dimensional (2D) vector of doubles containing the correlation values between the given 'Bacteria'
/// objects.
/// </returns>
vector<vector<double>> CompareAllBacteriaSequential()
{
	// Initialise a 2D vector to size of number of bacteria
	vector<vector<double>> corrVector(number_bacteria);

	// Instantiate array of pointers to 'Bacteria' objects
	Bacteria** b = new Bacteria * [number_bacteria];

	// Iterate through all bacteria and load them
	for (int i = 0; i < number_bacteria; i++)
	{
		b[i] = new Bacteria(bacteria_name[i]);
	}

	// Iterate through all 'Bacteria' objects
	for (int i = 0; i < number_bacteria - 1; i++) {
		// Initialise the row of columns to store correlation values
		corrVector[i] = vector<double>(number_bacteria, 0.0);

		// Calculate the correlation value and store it in the vector
		for (int j = i + 1; j < number_bacteria; j++) {
			double correlation = CompareBacteria(b[i], b[j]); // Calculate the correlation
			corrVector[i][j] = correlation; // Store in the vector
		}
	}

	// Return the resulting 2D vector containing the results
	return corrVector;
}

/// <summary>
///	Method to run timed tests on the parallel CompareAllBacteria() method.
/// </summary>
/// <param name="numTests">
/// The number of tests to run.
/// </param>
/// <param name="numThreads">
/// The number of threads to use.
/// </param>
void RunTests(int numTests, int numThreads) {
	SetThreads(numThreads); // Set the number of threads
	double* durations = new double[numTests]; // Create an array to store the test times

	for (int i = 0; i < numTests; i++) {
		printf("Starting test %d of %d.\n", i + 1, numTests);

		auto start = chrono::high_resolution_clock::now(); // Fetch the start time
		Init();
		ReadInputFile("list.txt");
		vector<vector<double>> resultVec = CompareAllBacteria();
		auto end = chrono::high_resolution_clock::now(); // Fetch the end time
		chrono::duration<double> elapsed = end - start; // Calculate the total time elapsed

		vector<double> parResult = MakeVectorOneDimension(resultVec);
		vector<double> seqResult = ReadVectorFromFile("./correlation.txt"); // Read sequential results from file

		// Verify the results are the same
		CompareResults(seqResult, parResult);

		// Clear the vectors
		resultVec.clear();
		parResult.clear();
		seqResult.clear();

		durations[i] = elapsed.count(); // Store the test's duration
	}

	// Print the test results to the terminal
	for (int i = 0; i < numTests; i++) {
		printf("%.10f\n", durations[i]);
	}

	delete[] durations; // Clear the array of durations
}

/// <summary>
///	Method to run timed tests on the sequential CompareAllBacteriaSequential() method.
/// </summary>
/// <param name="numTests">
/// The number of tests to run.
/// </param>
void RunSequentialTest(int numTests) {
	double* durations = new double[numTests];

	for (int i = 0; i < numTests; i++) {
		printf("Starting test %d of %d.\n", i + 1, numTests);
		auto start = chrono::high_resolution_clock::now(); // Fetch the start time

		Init();
		ReadInputFile("list.txt");
		vector<vector<double>> resultVec = CompareAllBacteriaSequential();

		auto end = chrono::high_resolution_clock::now(); // Fetch the end time
		chrono::duration<double> elapsed = end - start; // Calculate the total time elapsed
		resultVec.clear(); // Clear the vector

		durations[i] = elapsed.count(); // Store the tests duration
	}

	// Print the test results to the terminal
	for (int i = 0; i < numTests; i++) {
		printf("%.10f\n", durations[i]);
	}

	delete[] durations; // Clear the array of durations
}

/// <summary>
/// Main method for the running of the program.
/// </summary>
/// <param name="argc">
/// The number of arguments passed.
/// </param>
/// <param name="argv">
/// Array of pointers to arrays of character objects
/// </param>
/// <returns>
/// Implicit return.
/// </returns>
int main(int argc, char* argv[])
{
	RunTests(1, 8);
	//RunSequentialTest(10);

	return 0; // Exit program
}

/// ------------------------------------------------------------------------------------------------------------- ///
///												END OF FILE														  ///
/// ------------------------------------------------------------------------------------------------------------- ///