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
using namespace std;

/// -----------------------------------------------------------------------------------------------------------------
/// CVTree [TIDY] Version
/// -----------------------------------------------------------------------------------------------------------------
/*
int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

void Init()
{
	M2 = 1;
	for (int i = 0; i < LEN - 2; i++)	// M2 = AA_NUMBER ^ (LEN-2);
		M2 *= AA_NUMBER;
	M1 = M2 * AA_NUMBER;		// M1 = AA_NUMBER ^ (LEN-1);
	M = M1 * AA_NUMBER;			// M  = AA_NUMBER ^ (LEN);
}

class Bacteria
{
private:
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
		for (int i = 0; i < LEN - 1; i++) // Can parallelise
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
	long* vector;

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
		fclose(bacteria_file);
	}

	~Bacteria()
	{
		delete vector;
		delete second;
	}

	double stochastic_compute(long i)
	{
		double p1 = (double)second[i / AA_NUMBER] / (total + complement);
		double p2 = (double)one_l[i % AA_NUMBER] / total_l;
		double p3 = (double)second[i % M1] / (total + complement);
		double p4 = (double)one_l[i / M1] / total_l;
		return total * (p1 * p2 + p3 * p4) / 2;
	}
};

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

// Hotspot identified by VTune
double CompareBacteria(Bacteria* b1, Bacteria* b2)
{
	double correlation = 0;
	double vector_len1 = 0;
	double vector_len2 = 0;

	for (long i = 0; i < M; i++) // Split across threads
	{
		double stochastic1 = b1->stochastic_compute(i);
		double t1;
		if (stochastic1 > EPSILON)
			t1 = (b1->vector[i] - stochastic1) / stochastic1;
		else
			t1 = 0;
		vector_len1 += (t1 * t1);

		double stochastic2 = b2->stochastic_compute(i);
		double t2;
		if (stochastic2 > EPSILON)
			t2 = (b2->vector[i] - stochastic2) / stochastic2;
		else
			t2 = 0;
		vector_len2 += (t2 * t2);

		correlation = correlation + t1 * t2;
	}

	return correlation / (sqrt(vector_len1) * sqrt(vector_len2));
}

void CompareAllBacteria()
{
	// Create correlation array to store results
	double[] corr_arr = new double[number_bacteria - 1]

	for (int i = 0; i < number_bacteria - 1; i++)
	{
		Bacteria* b1 = new Bacteria(bacteria_name[i]);

		for (int j = i + 1; j < number_bacteria; j++) // Can convert this to linear (replace j with (i + 1))
		{
			Bacteria* b2 = new Bacteria(bacteria_name[j]);
			double correlation = CompareBacteria(b1, b2);
			printf("%03d %03d -> %.10lf\n", i, j, correlation); // Replace print with store bacteria
			delete b2;
		}
		delete b1;
	}
}

int main(int argc, char* argv[])
{
	time_t t1 = time(NULL);

	Init();
	ReadInputFile("list.txt");
	CompareAllBacteria();

	time_t t2 = time(NULL);
	printf("time elapsed: %lld seconds\n", t2 - t1);
	return 0;
}
*/
/// -----------------------------------------------------------------------------------------------------------------
/// CVTree [Improved] Version | Best Sequential Version
/// -----------------------------------------------------------------------------------------------------------------
int number_bacteria;
char** bacteria_name;
long M, M1, M2;
short code[27] = { 0, 2, 1, 2, 3, 4, 5, 6, 7, -1, 8, 9, 10, 11, -1, 12, 13, 14, 15, 16, 1, 17, 18, 5, 19, 3 };
#define encode(ch)		code[ch-'A']
#define LEN				6
#define AA_NUMBER		20
#define	EPSILON			1e-010

/// -----------------------------------------------------------------------------------------------------------------
/// New Methods
/// -----------------------------------------------------------------------------------------------------------------
/// <summary>
/// Method to write results to an output file.
/// </summary>
/// <param name="arr">
/// 
/// </param>
void WriteToFile(vector<double> arr) {
	std::ofstream outputFile("./correlation.txt"); // Set output file

	//
	for (const auto& value : arr) {
		outputFile << value << "\n";
	}
}

/// <summary>
/// Method to read results from an input file.
/// </summary>
/// <returns>
/// 
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
/// Method to compare the results of the sequential 
/// </summary>
/// <param name="seqCorrArr">
/// 
/// </param>
/// <param name="parCorrArr">
/// 
/// </param>
/// <returns>
/// 
/// </returns>
bool CompareResults(vector<double> seqCorrArr, vector<double> parCorrArr) {
	/*if (seqCorrArr == parCorrArr) {
		printf("Sequential and Parallel results match.");
		return true;
	}
	else {
		printf("Sequential and Parallel results do not match.");
		return false;
	}*/
	int count = 0;

	for (int i = 0; i < seqCorrArr.size(); i++) {
		//printf("Sequential: %f | Parallel: %f\n", seqCorrArr.at(i), parCorrArr.at(i));
		if (seqCorrArr.at(i) != parCorrArr.at(i)) {
			printf("Mismatch at %d: %f and %f.", i, seqCorrArr.at(i), parCorrArr.at(i));
			count++;
		}
	}

	/// TODO:
	///		[ ] Fix error where matching values are incorrectly listed as not matching
	///		[ ] Verify that the program will return true when comparing two matching values
	///		[ ] Extract program into a seperate project or even just file to be used later
	///				as a check & balance for further parallelisations.

	if (count > 0) {
		return false;
	}
	else {
		return true;
	}
}

/// -----------------------------------------------------------------------------------------------------------------
/// Unchanged Methods
/// -----------------------------------------------------------------------------------------------------------------
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

		for (long i = 0; i < M; i++) // Can be parallelised (LOW PRIO)
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
		for (long i = 0; i < M; i++) // Can be parallelised (LOW PRIO)
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

	for (long i = 0; i < number_bacteria; i++) // Can be parallelise (MEDIUM PRIO)
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

/// -----------------------------------------------------------------------------------------------------------------
/// Modified Methods
/// -----------------------------------------------------------------------------------------------------------------
/// <summary>
/// 
/// </summary>
vector<double> CompareAllBacteria()
{
	vector<double> result;
	int count = 0;
	Bacteria** b = new Bacteria * [number_bacteria];
	for (int i = 0; i < number_bacteria; i++) // Can be parallelised (HIGH PRIO)
	{
		//printf("load %d of %d\n", i + 1, number_bacteria); // Load information not important
		b[i] = new Bacteria(bacteria_name[i]);
	}

	for (int i = 0; i < number_bacteria - 1; i++) // Can be parallelised (HIGH PRIO)
	{
		for (int j = i + 1; j < number_bacteria; j++) // Can be converted to linear loop
		{
			//printf("%2d %2d -> ", i, j); // TODO: Remove this
			double correlation = CompareBacteria(b[i], b[j]);
			result.push_back(correlation);
			//printf("%.20lf\n", correlation); // TODO: Remove this
		}
	}
	printf("Array Size: %d\n", result.size());
	return result;
}

/// <summary>
/// 
/// </summary>
/// <param name="argc"></param>
/// <param name="argv"></param>
/// <returns></returns>
int main(int argc, char* argv[])
{
	time_t t1 = time(NULL);

	Init();
	ReadInputFile("list.txt");
	vector<double> result = CompareAllBacteria();
	WriteToFile(result);
	vector<double> checkResult = ReadFromFile();
	bool match = CompareResults(result, checkResult);
	if (match) {
		printf("The two vectors match.\n");
	}
	else {
		printf("The two vectors do not match.\n");
	}
	time_t t2 = time(NULL);
	printf("time elapsed: %lld seconds\n", t2 - t1);
	return 0;
}


/// -----------------------------------------------------------------------------------------------------------------
/// TODO:
///		[X] Add method to write output to file (binary)
///		[ ] Add method to verify new output against binary file
///		[ ] Parallelise high priority tasks
///		[ ] Convert triangular loops to linear where possible
///		[ ] Parallelise medium priority tasks
///		[ ] Parallelise low priority tasks (if possible)
///		[ ] Profile at each stage 
/// -----------------------------------------------------------------------------------------------------------------