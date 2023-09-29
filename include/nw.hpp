#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <fstream>
#include <filesystem>
#include <chrono>
// Neddleman-Wunsch
using namespace std;

#define inf numeric_limits<int>::min();
const string RESET_COLOR = "\033[0m";
const string RED_TEXT = "\033[31m";
const string GREEN_TEXT = "\033[32m";
const string YELLOW_TEXT = "\033[33m";
const string BLUE_TEXT = "\033[34m";
const string MAGENTA_TEXT = "\033[35m";
const string CYAN_TEXT = "\033[36m";
const string BOLD_TEXT = "\033[1m";
const int penaltyScore = 2;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<const string, string> pss;
typedef long long ll;
typedef unordered_map<string,string> MultipleAlignment;

inline int max3(int a, int b, int c)
{
	return max(max(a, b), c);
}

struct Alignment
{
	pair<string, string> s;
	pair<string, string> t;
};

void printAlignment(Alignment &a)
{
	cout << a.s.first << ": " << a.s.second << '\n';
	cout << a.t.first << ": " << a.t.second << '\n';
}

Alignment buildAlignment(vector<int> &path, pss &s, pss &t, int tam)
{
	int pathLength = tam - 1;
	int ti = 0;
	int si = 0;

	Alignment aligment;
	aligment.s.first = s.first;
	aligment.t.first = t.first;

	for (int dir = pathLength; dir >= 0; dir--)
	{
		if (path[dir] == 2)
			aligment.s.second += "-";
		else
			aligment.s.second += s.second[si++];
	}

	for (int dir = pathLength; dir >= 0; dir--)
	{
		if (path[dir] == 3)
			aligment.t.second += "-";

		else
			aligment.t.second += t.second[ti++];
	}
	return aligment;
}

void traceback(int i, int j, vvi &dp, pss &s, pss &t, vector<int> &path, ll &count, int &index, vector<Alignment> &aligments)
{
	if (i == 0 && j == 0)
	{
		Alignment aligment = buildAlignment(path, s, t, index);
		aligments.push_back(aligment);
		return;
	}
	else
	{
		int diag, up, left, max;
		bool match = s.second[i - 1] == t.second[j - 1];

		diag = up = left = inf;

		if (i > 0)
		{
			if (j > 0)
			{
				diag = dp[i - 1][j - 1] + (match ? +1 : -1);
				left = dp[i - 1][j] - penaltyScore;
				up = dp[i][j - 1] - penaltyScore;
			}
			else
			{
				left = dp[i - 1][j] - penaltyScore;
			}
		}
		else
		{
			up = dp[i][j - 1] - penaltyScore;
		}

		max = max3(diag, up, left);

		if (diag == max)
		{
			path[index++] = !match;
			traceback(i - 1, j - 1, dp, s, t, path, count, index, aligments);
			index--;
		}
		if (up == max)
		{
			path[index++] = 2;
			traceback(i, j - 1, dp, s, t, path, count, index, aligments);
			index--;
		}
		if (left == max)
		{
			path[index++] = 3;
			traceback(i - 1, j, dp, s, t, path, count, index, aligments);
			index--;
		}
	}
}

void show(vvi &dp)
{
	for (auto v : dp)
	{
		for (auto i : v)
		{
			cout << i << ' ';
		}
		cout << '\n';
	}
}

int nw(pss &s, pss &t, vector<vector<int>> &dp)
{
	int rows, cols;
	rows = s.second.size();
	cols = t.second.size();

	dp.clear();
	dp.resize(rows + 1, vector<int>(cols + 1, 0));

	for (int i = 1; i <= rows; i++)
		dp[i][0] = -penaltyScore * i;

	for (int j = 1; j <= cols; j++)
		dp[0][j] = -penaltyScore * j;

	for (int i = 1; i <= rows; i++)
		for (int j = 1; j <= cols; j++)
			dp[i][j] = max3(
					dp[i - 1][j - 1] + (s.second[i - 1] == t.second[j - 1] ? +1 : -1),
					dp[i - 1][j] - penaltyScore,
					dp[i][j - 1] - penaltyScore);

	return dp[rows][cols];
}

int nw(pair<string, string> &s, pss &t, vector<vector<int>> &dp)
{
	pss sc(s.first, s.second);
	return nw(sc, t, dp);
}

vector<Alignment> getAlignments(pss &s, pss &t, vector<vector<int>> &dp)
{
	int rows, cols;
	vector<Alignment> aligments;

	rows = s.second.size();
	cols = t.second.size();

	vector<int> tmpPath(rows + cols);
	ll count = 0;
	int index = 0;

	traceback(rows, cols, dp, s, t, tmpPath, count, index, aligments);

	return aligments;
}
vector<Alignment> getAlignments(pair<string, string> &s, pss &t, vector<vector<int>> &dp)
{
	pss sc(s.first, s.second);
	return getAlignments(sc, t, dp);
}

void printAlignments(vector<Alignment> &aligments, pss &s, pss &t, vvi &dp)
{
	int rows, cols;
	rows = s.second.size();
	cols = t.second.size();

	cout << BOLD_TEXT << GREEN_TEXT << "[ INFO | Original Sequences ]\n"
			 << RESET_COLOR;

	cout << s.first << ": " << s.second << '\n';
	cout << t.first << ": " << t.second << "\n\n";

	cout << BOLD_TEXT << GREEN_TEXT << "[ INFO | Possible Alignments ]\n"
			 << RESET_COLOR;

	int c = 0;
	for (auto aligment : aligments)
	{
		cout << "[" << ++c << "]\n";
		printAlignment(aligment);
	}

	cout << BOLD_TEXT << GREEN_TEXT << "[ INFO | Summary]\n";
	cout << "max score: " << dp[rows][cols] << '\n';
	cout << "n of alignments: " << aligments.size() << RESET_COLOR << '\n';
}

void printAlignments(vector<Alignment> &aligments, pair<string, string> &s, pss &t, vvi &dp)
{
	pss sc(s.first, s.second);
	printAlignments(aligments, sc, t, dp);
}

vector<string> getFilenames(string &dir)
{
	vector<string> filenames;
	for (const auto &entry : filesystem::directory_iterator(dir))
	{
		if (entry.is_regular_file())
			filenames.push_back(entry.path());
	}
	return filenames;
}

unordered_map<string, string> getSequences(vector<string> &filenames)
{
	unordered_map<string, string> sequences;

	for (auto filename : filenames)
	{
		fstream file(filename);
		string word, sequence, name;

		file >> name;

		while (file >> word)
		{
			if (word.size() == 10)
				sequence += word;
		}
		sequences[name] = sequence;
	}
	return sequences;
}

void printSequences(unordered_map<string, string> &sequences)
{
	for (auto seq : sequences)
	{
		cout << "Name: " << seq.first << '\n';
		cout << "Length: " << seq.second.size() << '\n';
		cout << "Sequence: " << seq.second << "\n\n";
	}
}

// int main()
// {
// 	string dir = "../sequences";
// 	vector<string> filenames = getFilenames(dir);
// 	unordered_map<string, string> sequences = getSequences(filenames);
// 	printSequences(sequences);

// 	string s, t;
// 	s = sequences["Bacteria"];
// 	t = sequences["Sars-Cov"];
// 	// s = "AAAC";
// 	// t = "AGC";
// 	nw(s, t);
// }