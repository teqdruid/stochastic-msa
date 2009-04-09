#include "storage.h"

void parse(string& identifier, GeneticSymbols*& data,
	   size_t& size, istream& inp) {
    char buffer[128];
    size_t cap = 100;
    size = 0;
    data = (GeneticSymbols*)malloc(sizeof(GeneticSymbols) * cap);

    if (inp.peek() == '>') {
	char buffer[1024];
	inp.getline(buffer, 1024);
	identifier = string(&buffer[1]);
    }

    while (!inp.eof()) {
	inp.get(buffer, 128, '>');
	size_t bfill = inp.gcount();

	for (size_t i=0; i<bfill; i++) {
	    char c = buffer[i];
	    if (c == ' ' || c == '\t' || c == '\n' || c == '\r'
		|| c == 'N' || c == 'n')
		continue; 

	    if (size > cap) {
		cap *= 2;
		data = (GeneticSymbols*)realloc(data, sizeof(GeneticSymbols) * cap);
	    }

	    data[size++] = fromChar(c);
	}

	if (inp.eof() || inp.peek() == '>')
	    break;
    }

    data = (GeneticSymbols*)realloc(data, sizeof(GeneticSymbols) * size);
}


