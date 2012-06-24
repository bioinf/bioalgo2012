#ifndef _FASTA_READER_H_
#define _FASTA_READER_H_

#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream> 

class FASTAReader
{
public:
	bool IsOk() const
	{
		return inStream_ != 0;
	}

	explicit FASTAReader(const std::string & fileName): inStream_(fileName.c_str())
	{

	}

	void GetSequence(std::string & buf)
	{
		buf.clear();
		std::string temp;
		inStream_ >> temp;
		while(inStream_ >> temp)
		{
			buf += temp;
		}
	}
		
private:
	std::ifstream inStream_;
};


#endif 
