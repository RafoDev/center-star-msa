#include <iostream>
#include <algorithm>
#include "include/nw.hpp"
using namespace std;

void adjustAlignments(vector<pair<string, string>> &msa, vi adjustmentPositions)
{
  for (auto position : adjustmentPositions)
    for (auto &a : msa)
      a.second.insert(position, 1, '-');
}

void fillGaps(vector<pair<string, string>> &msa)
{
  int maxSize = 0;
  for (auto &a : msa)
    maxSize = max(int(a.second.size()), maxSize);

  for (auto &a : msa)
  {
    if (a.second.size() < maxSize)
      a.second.insert(a.second.end(), maxSize - a.second.size(), '-');
  }
}

void pairwiseAlignmentUtil(vector<pair<string, string>> &msa, Alignment &a, pair<string, string> &center)
{
  string *originalCenter = &msa[0].second;
  string *currCenter = &a.s.second;

  int originalCenterSize = originalCenter->size();
  int currCenterSize = currCenter->size();

  vi adjustmentPositions;

  int originalCenterGaps = 0;
  int currCenterGaps = 0;

  if (originalCenterSize != currCenterSize)
  {
    int i = 0, j = 0;
    while (i < originalCenterSize && j < currCenterSize)
    {
      if (originalCenter->at(i) != currCenter->at(j))
      {
        if (originalCenter->at(i) == '-')
        {
          i++;
        }
        else
        {
          adjustmentPositions.push_back(i);
          j++;
        }
      }
      else
      {
        i++;
        j++;
      }
    }
    adjustAlignments(msa, adjustmentPositions);
  }
  msa.push_back(a.t);
  fillGaps(msa);
}

void pairwiseAlignment(vector<Alignment> multipleAlignments, pair<string, string> &center)
{

  vector<pair<string, string>> msa;
  if (!multipleAlignments.empty())
  {
    cout << BOLD_TEXT << GREEN_TEXT << "[INFO]" << RESET_COLOR << " Alignments with Sc\n\n";

    for (auto a : multipleAlignments)
    {
      printAlignment(a);
      cout << '\n';
    }

    msa.push_back(multipleAlignments[0].s);
    msa.push_back(multipleAlignments[0].t);

    int multipleAlignmentsSize = multipleAlignments.size();

    cout << BOLD_TEXT << GREEN_TEXT << "[INFO]" << RESET_COLOR << " Pairwise Alignment\n";

    for (int i = 1; i < multipleAlignmentsSize; i++)
    {
      cout << "\n-- step " << i << ":\n\n";
      for (auto alg : msa)
        cout << alg.first << ": " << alg.second << '\n';

      cout << '\n';

      cout << multipleAlignments[i].s.first << ": " << multipleAlignments[i].s.second << '\n';
      cout << multipleAlignments[i].t.first << ": " << multipleAlignments[i].t.second << '\n';

      cout << "\n-- result: \n\n";

      pairwiseAlignmentUtil(msa, multipleAlignments[i], center);

      for (auto alg : msa)
        cout << alg.first << ": " << alg.second << '\n';
    }
    cout << '\n';
  }
}

void starMSA(unordered_map<string, string> &sequences)
{

  vector<vector<int>> dp;

  int nSequences = sequences.size();
  vector<vector<int>> scores(nSequences, vector<int>(nSequences, 0));

  int maxScore, i, j;
  i = j = 0;
  maxScore = inf;
  pair<string, string> center;

  for (auto &seqI : sequences)
  {
    int rowScores = 0;
    j = 0;
    for (auto &seqJ : sequences)
    {
      if (seqI.first != seqJ.first)
        scores[i][j] = nw(seqI, seqJ, dp);

      rowScores += scores[i][j];
      j++;
    }
    if (rowScores > maxScore)
    {
      maxScore = rowScores;
      center = seqI;
    }
    i++;
  }

  cout << BOLD_TEXT << GREEN_TEXT << "[INFO]" << RESET_COLOR << " Center\n\n";

  cout << center.first << ": " << center.second << "\n\n";

  vector<Alignment> multipleAlignments;
  unordered_map<string, string> msa;

  for (auto seqI : sequences)
  {
    if (seqI.first != center.first)
    {
      nw(center, seqI, dp);
      vector<Alignment> alignments = getAlignments(center, seqI, dp);
      // printAlignments(alignments, center, seqI, dp);
      multipleAlignments.push_back(alignments.back());
    }
  }

  pairwiseAlignment(multipleAlignments, center);
}

int main()
{
  unordered_map<string, string> testSequences;

  testSequences["s1"] = "ATTGCCATT";
  testSequences["s2"] = "ATGGCCATT";
  testSequences["s3"] = "ATCCAATTTT";
  testSequences["s4"] = "ATCTTCTT";
  testSequences["s5"] = "ACTGACC";

  unordered_map<string, string> BRCASequencesF;

  BRCASequencesF["s1-F"] = "TGCCGGCAGGGATGTGCTTG";
  BRCASequencesF["s2-F"] = "GTTTAGGTTTTTGCTTATGCAGCATCCA";
  BRCASequencesF["s3-F"] = "GGAAAAGCACAGAACTGGCCAACA";
  BRCASequencesF["s4-F"] = "GCCAGTTGGTTGATTTCCACCTCCA";
  BRCASequencesF["s5-F"] = "ACCCCCGACATGCAGAAGCTG";
  BRCASequencesF["s6-F"] = "TGACGTGTCTGCTCCACTTCCA";

  unordered_map<string, string> BRCASequencesR;

  BRCASequencesR["s1-R"] = "TGCTTGCAGTTTGCTTTCACTGATGGA";
  BRCASequencesR["s2-R"] = "TCAGGTACCCTGACCTTCTCTGAAC";
  BRCASequencesR["s3-R"] = "GTGGGTTGTAAAGGTCCCAAATGGT";
  BRCASequencesR["s4-R"] = "TGCCTTGGGTCCCTCTGACTGG";
  BRCASequencesR["s5-R"] = "GTGGTGCATTGATGGAAGGAAGCA";
  BRCASequencesR["s6-R"] = "AGTGAGAGGAGCTCCCAGGGC";

  unordered_map<string, string> BRCASequencesFR;

  BRCASequencesFR["s1"] = BRCASequencesF["s1-F"] + BRCASequencesR["s1-R"];
  BRCASequencesFR["s2"] = BRCASequencesF["s2-F"] + BRCASequencesR["s2-R"];
  BRCASequencesFR["s3"] = BRCASequencesF["s3-F"] + BRCASequencesR["s3-R"];
  BRCASequencesFR["s4"] = BRCASequencesF["s4-F"] + BRCASequencesR["s4-R"];
  BRCASequencesFR["s5"] = BRCASequencesF["s5-F"] + BRCASequencesR["s5-R"];
  BRCASequencesFR["s6"] = BRCASequencesF["s6-F"] + BRCASequencesR["s6-R"];

  starMSA(testSequences);
  starMSA(BRCASequencesF);
  starMSA(BRCASequencesR);
  starMSA(BRCASequencesFR);
}