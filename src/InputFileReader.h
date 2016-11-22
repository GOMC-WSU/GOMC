#pragma once
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

class InputFileReader
{
private:
	fstream fs;
	vector<string> & split(const string &s, char delim, vector<string> &elems);

public:
	bool readNextLine(vector<string> & str);
	void Open(string fileName);
	InputFileReader(string fileName);
	InputFileReader(void);
	~InputFileReader();
};