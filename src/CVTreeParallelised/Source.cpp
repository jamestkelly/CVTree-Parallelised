#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <vector>
#include <fstream>
#include <iterator>
#include <string>
#include <iostream>
#include <omp.h>
using namespace std;

/// ------------------------------------------------------------------------------------------------------------- ///
///									CVTree [Improved] Version | Parallelised									  ///
/// ------------------------------------------------------------------------------------------------------------- ///
int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010
#define THREAD_NUM		4

/// ------------------------------------------------------------------------------------------------------------- ///
///													New Methods													  ///
/// ------------------------------------------------------------------------------------------------------------- ///
/// <summary>
/// Method to write results to an output file. This output comes in the expected format of a ".txt" file with the
/// calculated results of the correlation values between any two given bacteria in the data set.
/// </summary>
/// <param name="arr">
/// A vector of doubles, containing the correlation values as calculated to be written to the output file.
/// </param>
void WriteToFile(vector<double> arr) {
	std::ofstream outputFile("./correlation.txt"); // Set output file

	// Iterate over all values in the array and write to file
	for (const auto& value : arr) {
		outputFile << value << "\n";
	}
}

/// <summary>
/// Method to read results from an input file. This input comes in the expected format of a ".txt" file with the
/// calculated results of the correlation values between any two given bacteria in the data set.
/// </summary>
/// <returns>
/// A vector of doubles, containing the correlation values as calculated to be read from the input file.
/// </returns>
vector<double> ReadFromFile() {
	double num = 0.0; // Initialise value to store line
	vector<double> inputResult; // Initialise vector to store results
	std::ifstream inputFile("./correlation.txt", std::ios::in); // Set input file

	// Check to see the file was opened correctly
	if (!inputFile.is_open()) {
		std::cerr << "There was a problem opening the input file.\n"; // Print error message
		exit(1); // Exit
	}

	// Store values from input file until end of file (EOF)
	while (inputFile >> num) {
		inputResult.push_back(num);
	}

	return inputResult; // Return the resulting vector
}

/// <summary>
/// Method to perform an element-wise comparison of two vectors of doubles. These vectors contain the
/// correlation values calculated from the comparisons of bacteria.
/// </summary>
/// <param name="seqCorrArr">
/// The original correlation values as calculated previously through the use of the original best
/// sequential version of the program.
/// </param>
/// <param name="parCorrArr">
/// The correlation values as calculated in the current run of the program.
/// </param>
/// <returns>
/// Returns a true or false value depending on the whether the results match, element per element. True
/// is returned if all elements match to their corresponding values in the other vector, and False is
/// returned otherwise.
/// </returns>
bool CompareResults(vector<double> seqCorrArr, vector<double> parCorrArr) {
	double epsilon = 1e-5; // Set epsilon for maximum difference
	double difference = 0.0; // Initialise difference variable to 0

	// If the two vectors are not the same length, return false
	if (parCorrArr.size() != seqCorrArr.size()) {
		return false;
	}

	// Iterate through all results
	for (int i = 0; i < seqCorrArr.size(); i++) {
		// Calculate the difference between the two calculated correlations
		difference = abs(seqCorrArr.at(i) - parCorrArr.at(i));

		// If there is a significant difference in correlation
		if (difference >= epsilon) {
			// Print where the mismatch occurred
			printf("Mismatch at %d: %.11f and %.11f.\n", i, seqCorrArr.at(i), parCorrArr.at(i));

			return false; // Break the loop and return false
		}
	}

	return true; // Default to return true for no mismatches found
}

/// <summary>
/// Method to print a message to the terminal, displaying whether or not the correlation values calculated
/// match element for element.
/// </summary>
/// <param name="resultMatch">
/// A true or false value corresponding to the accuracy of the parallel results calculated.
/// </param>
void ExitCheck(bool resultMatch) {
	if (resultMatch) {
		printf("The parallel results match the sequential correlation values.\nExiting program...");
	}
	else {
		printf("The parallel results do not match the sequential correlation values.\nExiting program...");
	}
}

/// ------------------------------------------------------------------------------------------------------------- ///
///													Unchanged Methods											  ///
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
	/// <summary>
	/// 
	/// </summary>
	private:
		long* vector;
		long* second;
		long one_l[AA_NUMBER];
		long indexs;
		long total;
		long total_l;
		long complement;

	/// <summary>
	/// 
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
	/// 
	/// </summary>
	/// <param name="buffer"></param>
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
	/// 
	/// </summary>
	/// <param name="ch"></param>
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

	/// <summary>
	/// 
	/// </summary>
	public:
		long count;
		double* tv;
		long* ti;

	/// <summary>
	/// 
	/// </summary>
	/// <param name="filename"></param>
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

		for (long i = 0; i < M; i++) // Can be parallelised (MEDIUM PRIO)
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
		for (long i = 0; i < M; i++) // Can be parallelised (MEDIUM PRIO)
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
/// <param name="input_name">
/// 
/// </param>
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

	for (long i = 0; i < number_bacteria; i++) // Can be parallelised MAYBE (LOW PRIO)
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
/// <param name="b1">
/// 
/// </param>
/// <param name="b2">
/// 
/// </param>
/// <returns>
/// 
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

/// ------------------------------------------------------------------------------------------------------------- ///
///													Modified Methods											  ///
/// ------------------------------------------------------------------------------------------------------------- ///
/// <summary>
/// TODO:
///		[ ] Add description of method and an explanation of the changes made to the source code.
/// </summary>
vector<double> CompareAllBacteria()
{
	// Initialise variables
	vector<double> result; // Vector to store correlation values
	int count = 0; // Initialise outer-loop counter
	int num_iter = number_bacteria - 1; // Set number of iterations to the number of bacteria
	Bacteria** b = new Bacteria * [number_bacteria];

	// Load the bacteria objects
	for (int i = 0; i < number_bacteria; i++) // Can be parallelised (HIGH PRIO)
	{
		b[i] = new Bacteria(bacteria_name[i]);
	}

	// Iterate through all bacteria and calculate correlation values.
	for (int i = 0; i < num_iter; i++) {
		result.push_back(CompareBacteria(b[count], b[i + 1]));

		// If the loop is about to reach the number of iterations
		if (i == num_iter - 1) {
			i = count; // Set 'i' to the counter's value
			count++; // Increment the counter
		}
	}

	return result; // Return the resulting vector of correlation values
}

/// <summary>
/// TODO:
///		[ ] Add description of method and an explanation of the changes made to the source code.
/// </summary>
/// <param name="argc">
/// 
/// </param>
/// <param name="argv">
/// 
/// </param>
/// <returns>
///		[ ] Add description of what the method returns.
/// </returns>
int main(int argc, char* argv[])
{
	time_t t1 = time(NULL);

	Init();
	ReadInputFile("list.txt");
	vector<double> parResult = CompareAllBacteria();
	time_t t2 = time(NULL);
	printf("time elapsed: %lld seconds\n", t2 - t1);


	// WriteToFile(result); // Write the results to file for verification
	vector<double> seqResult = ReadFromFile(); // Read sequential results from file

	// Check the results match
	ExitCheck(CompareResults(seqResult, parResult));
	
	return 0; // Exit
}

/// ------------------------------------------------------------------------------------------------------------- ///
///														EOF														  ///
/// ------------------------------------------------------------------------------------------------------------- ///