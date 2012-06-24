#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <iostream>
#include <algorithm>
#include <ANN/ANN.h>

//Library for KNN can be found at: http://www.cs.umd.edu/~mount/ANN/

#include "fasta.h"

const std::string abc = "ACGT";
const std::string rcomp = "TGCA";

char FindRComp(char ch)
{
	return rcomp[std::find(abc.begin(), abc.end(), ch) - abc.begin()];
}

struct Record
{
	int start;
	int end;
	int classId;	
	std::vector<float> coord;

	Record() {}
	Record(int start, int end, int classId, std::vector<float> coord): 
		start(start), end(end), classId(classId), coord(coord) {}
};

const int CHR_COUNT = 22;

template<class It>
	void GenerateChromosomeList(int n, It out)
	{
		for(int i = 1; i <= n; i++)
		{
			std::stringstream ss;
			ss << "chr" << i;
			*out++ = ss.str();
		}
		
		*out++ = "chrX";
		*out++ = "chrY";
	}

template<class It>
	void GoKMer(int now, int k, const std::string & abc, std::string & buf, It out)
	{
		if(now == k)
		{
			*out++ = buf;
		}
		else
		{
			for(size_t i = 0; i < abc.size(); i++)
			{
				buf[now] = abc[i];
				GoKMer(now + 1, k, abc, buf, out);
			}
		}
	}

template<class It>
	void GenerateKmerList(int k, const std::string & abc, It out)
	{
		std::string buf(k, ' ');
		GoKMer(0, k, abc, buf, out);
	}


void MakeRCompTranslationTable(const std::vector<std::string> & kmer, std::vector<size_t> & ret)
{
	ret.resize(kmer.size());
	std::string buf(kmer[0].size(), ' ');
	for(size_t i = 0; i < kmer.size(); i++)
	{
		std::transform(kmer[i].rbegin(), kmer[i].rend(), buf.begin(), FindRComp);
		ret[i] = find(kmer.begin(), kmer.end(), buf) - kmer.begin();
	}
}

int GetIndex(const std::string & str, const std::vector<std::string> & kmer)
{
	return static_cast<int>(std::lower_bound(kmer.begin(), kmer.end(), str) - kmer.begin());
}

void MapToVector(const std::string & str, int k, int size, 
	const std::vector<std::string> & kmer, std::vector<float> & ret)
{
	ret.assign(size, 0);
	for(size_t i = 0; i < str.size() - k + 1; i++)
	{
		std::string temp(str.begin() + i, str.begin() + i + k);
		if(temp.find('N') == temp.npos)
		{
			int index = GetIndex(temp, kmer);
			ret[index] += 1.0 / (str.size() - k);
		}
	}
}

template<class It>
	void ReadObjects(int k, bool train, std::vector<std::string> & kmer, const std::string & abc,
		const std::string & genomeFile, const std::string & classFile, It out, int limit = -1)
	{
		std::string genome;
		FASTAReader reader(genomeFile);
		if(!reader.IsOk())
		{
			throw std::string("Cant open file with the genome: ") + genomeFile;
		}

		std::ifstream classIn(classFile.c_str());
		if(!classIn)
		{
			throw std::string("Cant open file with the classes: ") + classFile;
		}
				
		int end;
		int start;
		int classId = -1;
		int size = 1;
		for(int i = 0; i < k; i++)
		{
			size *= static_cast<int>(abc.size());
		}

		std::vector<float> temp;
		reader.GetSequence(genome);		
		for(int count = 0; ; count++)
		{
			classIn >> start >> end;
			if(train)
			{
				classIn >> classId;
			}

			if(!classIn)
			{
				break;
			}

			std::string str(genome.begin() + start - 1, genome.begin() + end - 1);
			MapToVector(str, k, size, kmer, temp);
			*out++ = Record(start, end, classId, temp);

			if(limit != -1 && count > limit)
			{
				break;
			}
		}
	}

int Query(const std::vector<Record> & trainSet, int nn, ANNkd_tree * kdTree, Record & q)
{
	static std::vector<int> nnIndex(nn);
	static std::vector<double> nnDists(nn);

	std::map<int, int> cnt;
	std::vector<double> temp(q.coord.begin(), q.coord.end());
	kdTree->annkSearch(&temp[0], nn, &nnIndex[0], &nnDists[0], 0.00);
	for(int i = 0; i < nn; i++)
	{
		int classId = trainSet[nnIndex[i]].classId;
		if(cnt.find(nnIndex[i]) == cnt.end())
		{
			cnt[classId] = 0;
		}

		cnt[classId] += 1;
	}

	int maxScore = 0;
	for(std::map<int, int>::iterator it = cnt.begin(); it != cnt.end(); it++)
	{
		maxScore = std::max(maxScore, it->second);
	}

	for(std::map<int, int>::iterator it = cnt.begin(); it != cnt.end(); it++)
	{
		if(it->second == maxScore)
		{
			return it->first;
		}
	}

	return -1;
}

std::pair<int, double> Query1(const std::vector<Record> & trainSet, ANNkd_tree * kdTree, Record & q)
{
	int nnIndex;
	double nnDists;

	std::vector<double> temp(q.coord.begin(), q.coord.end());
	kdTree->annkSearch(&temp[0], 1, &nnIndex, &nnDists, 0.00);
	return std::make_pair(trainSet[nnIndex].classId, nnDists);
}

ANNkd_tree * BuildTree(const std::vector<Record> & record)
{
	ANNpointArray dataPoints = annAllocPts(record.size(), record[0].coord.size());
	for(size_t i = 0; i < record.size(); i++)
	{
		std::copy(record[i].coord.begin(), record[i].coord.end(), dataPoints[i]);
	}

	return new ANNkd_tree(dataPoints, record.size(), record[0].coord.size());
}

void Translate(std::vector<float> & src, std::vector<size_t> ttable, std::vector<float> & ret)
{
	ret.resize(src.size());
	for(size_t i = 0; i < ret.size(); i++)
	{
		ret[ttable[i]] = src[i];
	}
}

int main(int argc, char * argv[])
{
	if(argc == 7)
	{				
		std::vector<Record> testSet;
		std::vector<Record> trainSet;

		int k = atoi(argv[1]);
		int nn = atoi(argv[2]);
		std::vector<std::string> kmer;
		std::string genomeDir = argv[3];
		std::string trainDir = argv[4];
		std::string testDir = argv[5];
		std::string ansDir = argv[6];
		std::vector<std::string> chromosome;
		GenerateKmerList(k, abc, std::back_inserter(kmer));
		GenerateChromosomeList(CHR_COUNT, std::back_inserter(chromosome));
		std::vector<size_t> translationTable;
		MakeRCompTranslationTable(kmer, translationTable);			
		Record rcomp(0, 0, 0, std::vector<float>(kmer.size()));
		try
		{
			std::cout << "Loading train data" << std::endl;
			for(size_t i = 0; i < chromosome.size(); i++)
			{
				std::cout << "Loading chr No. " << i << std::endl;
				std::string classFile = trainDir + "/" + chromosome[i];
				std::string genomeFile = genomeDir + "/" + chromosome[i] + ".fasta";
				ReadObjects(k, true, kmer, abc, genomeFile, classFile, std::back_inserter(trainSet));
			}

			ANNkd_tree * tree = BuildTree(trainSet);
			std::cout << "Starting classification" << std::endl;

			int total = 0;
			for(size_t i = 0; i < chromosome.size(); i++)
			{
				std::cout << "Loading chr No. " << i << std::endl;
				std::string ansFile = ansDir + "/" + chromosome[i];
				std::string classFile = testDir + "/" + chromosome[i];
				std::string genomeFile = genomeDir + "/" + chromosome[i] + ".fasta";
				std::ofstream ansOut(ansFile.c_str());	

				testSet.clear();
				ReadObjects(k, false, kmer, abc, genomeFile, classFile, std::back_inserter(testSet));
				for(size_t j = 0; j < testSet.size(); j++)
				{	
					Translate(testSet[j].coord, translationTable, rcomp.coord);
					std::pair<int, double> q1 = Query1(trainSet, tree, testSet[j]);
					std::pair<int, double> q2 = Query1(trainSet, tree, rcomp);
					ansOut << testSet[j].start << ' ' << testSet[j].end << ' ' << (q1.second < q2.second ? q1.first : q2.first) << std::endl;
					if(++total % 10000 == 0)
					{
						std::cout << total << std::endl;
					}
				}
			}

			delete tree;
		}
		catch(std::string & buf)
		{
			std::cerr << buf;
			return 1;
		}
	}
	else
	{
		std::cout << "Usage: Classify <k-mer> <nn> <genome dir> <train dir> <test dir> <answer dir>" << std::endl;
	}

	return 0;
}