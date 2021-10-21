#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace std;

/// ------------------------------------------------------------------------------------------------------------- ///
///											CVTree [Improved] | Parallelised									  ///
/// ------------------------------------------------------------------------------------------------------------- ///
/// Author: Jim Tran Kelly | N9763686																			  ///
///																												  ///
/// Description: This script has been updated and parallelised for the purposes of CAB401 High Performance and    ///
/// Parallel Computing at the Queensland University of Technology (QUT), Semester 2, 2021. The script in its      ///
/// original form is a program to take a file as input through the command-line, containing the number of strings ///
/// in the first line and following that several bacteria names. These names are then concatenated with .faa	  ///
/// matching to files in the /data/ directory with gene names and their respective gene represented as DNA. The   ///
/// program collects the "kmer" or subset of an entire gene and then compares that against other bacteria "kmers" ///
/// and finds the correlation between any given bacteria. This application was found on the Blackboard Assignment ///
/// tab for CAB401 as one of the provided projects that is available for parallelization.						  ///
///																												  ///
/// Version: 3.1																								  ///
/// ------------------------------------------------------------------------------------------------------------- ///
int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

/// ------------------------------------------------------------------------------------------------------------- ///
///													New Methods													  ///
/// ------------------------------------------------------------------------------------------------------------- ///

/// <summary>
/// 
/// </summary>
/// <param name="corrArr">
/// 
/// </param>
/// <param name="fileName">
/// 
/// </param>
void WriteToFile(vector<vector<double>> corrArr, string fileName) {
	std::ofstream outputFile(fileName); // Set output file
	
	for (int i = 0; i < corrArr.size(); i++) {
		for (int j = 0; j < corrArr.at(i).size(); j++) {
			outputFile << std::defaultfloat 
				<< std::setprecision(std::numeric_limits<long double>::digits10) 
				<< corrArr[i][j] << "\n";
		}
	}
}

/// <summary>
/// 
/// </summary>
/// <param name="fileName">
/// 
/// </param>
/// <returns>
/// 
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
/// 
/// </summary>
/// <param name="corrVector">
/// 
/// </param>
/// <returns>
/// 
/// </returns>
vector<double> MakeVectorOneDimension(vector<vector<double>> corrVector) {
	vector<double> resultVector; // Initialise one dimensional vector of doubles

	// Iterate through all values in the input vector
	for (int i = 0; i < corrVector.size(); i++) {
		for (int j = 0; j < corrVector[i].size(); j++) {
			resultVector.push_back(corrVector[i][j]); // Append to the output vector
		}
	}

	return resultVector; // Return the resulting one dimensional vector of doubles
}

/// <summary>
/// 
/// </summary>
/// <param name="seqResult">
/// 
/// </param>
/// <param name="parResult">
/// 
/// </param>
void CompareResults(vector<double> seqResult, vector<double> parResult) {
	double difference = 0.0; // Initialise difference variable to 0
	int mismatchCount = 0; // Initialise a counter for the number of mismatches

	// Iterate through all pairs of values
	for (int i = 0; i < seqResult.size(); i++) {
		difference = abs(seqResult[i] - parResult[i]); // Calculate the difference

		// If there is a significant difference in correlation
		if (difference >= EPSILON) {
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

/// ------------------------------------------------------------------------------------------------------------- ///
///												Unchanged Methods												  ///
/// ------------------------------------------------------------------------------------------------------------- ///

/// <summary>
/// 
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
/// 
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
	double* tv;
	long* ti;

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
		double* t = new double[M];

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
				t[i] = (vector[i] - stochastic) / stochastic;
				count++;
			}
			else
				t[i] = 0;
		}

		delete second_div_total;
		delete vector;
		delete second;

		tv = new double[count];
		ti = new long[count];

		int pos = 0;
		for (long i = 0; i < M; i++)
		{
			if (t[i] != 0)
			{
				tv[pos] = t[i];
				ti[pos] = i;
				pos++;
			}
		}
		delete t;

		fclose(bacteria_file);
	}
};

/// <summary>
/// 
/// </summary>
/// <param name="input_name"></param>
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
/// 
/// </summary>
/// <param name="b1"></param>
/// <param name="b2"></param>
/// <returns></returns>
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

/// ------------------------------------------------------------------------------------------------------------- ///
///												Modified Methods												  ///
/// ------------------------------------------------------------------------------------------------------------- ///

/// <summary>
/// 
/// </summary>
/// <returns>
/// 
/// </returns>
vector<vector<double>> CompareAllBacteria()
{
	// Initialise counter variables
	int rowCount = 0;
	int colCount = 1;
	int total_iter = 0;

	// Initialise a 2D vector to size of number of bacteria
	vector<vector<double>> corrVector(number_bacteria);

	Bacteria** b = new Bacteria * [number_bacteria];

	// Iterate through all bacteria and load them
	for (int i = 0; i < number_bacteria; i++)
	{
		//printf("load %d of %d\n", i + 1, number_bacteria);
		b[i] = new Bacteria(bacteria_name[i]);
	}

	// Iterate to create the 2D vector to store correlation values
	for (int i = 0; i < number_bacteria - 1; i++) {
		// Initialise the 'row' of values
		corrVector[i] = vector<double>(number_bacteria, 0.0);

		for (int j = 0; j < number_bacteria; j++) {
			total_iter++; // Increment total iteration counter
		}
	}

	// Iterate through all possible comparisons
	for (int i = 0; i < total_iter / 2; i++) {
		//printf("Iteration: %d | %2d %2d -> ", i, rowCount, colCount); // TODO: Remove this after debugging
		//double correlation = CompareBacteria(b[rowCount], b[colCount]);
		//printf("%.20lf\n", correlation); // TODO: Remove this after debugging

		corrVector[rowCount][colCount] = CompareBacteria(b[rowCount], b[colCount]);

		// If the 'column' count reaches the number of bacteria
		if (colCount == number_bacteria - 1) {
			rowCount++; // Increment the 'row' count
			colCount = rowCount + 1; // Set the 'column' count to the 'row' count plus one
		}
		else {
			colCount++; // Increment the 'column' count
		}
	}

	// Return the resulting 2D vector containing the results
	return corrVector;
}

/// <summary>
/// 
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(int argc, char* argv[])
{
	printf("Starting operations...\n");
	time_t t1 = time(NULL);

	Init();
	ReadInputFile("list.txt");
	vector<vector<double>> resultVec = CompareAllBacteria();
	time_t t2 = time(NULL);
	printf("time elapsed: %lld seconds\n", t2 - t1);
	printf("Correlation operations complete...\n");

	//WriteToFile(resultVec, "./correlation.txt"); // Write the sequential results to output file
	printf("Comparing results...\n");
	vector<double> seqResult = ReadVectorFromFile("./correlation.txt"); // Read sequential results from file
	vector<double> parResult = MakeVectorOneDimension(resultVec); // Convert the vector to one dimension
	
	// Verify the results are the same
	CompareResults(seqResult, parResult); 

	return 0;
}