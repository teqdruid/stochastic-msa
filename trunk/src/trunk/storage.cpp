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
	    if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
		continue; 

	    if (size >= cap) {
		cap *= 2;
		data = (GeneticSymbols*)realloc(data, sizeof(GeneticSymbols) * cap);
	    }

	    try {
		data[size++] = fromChar(c);
	    } catch (MsaException* e) {
		cerr << "Error reading character: " << c << ". Skipping sequence "
		     << identifier << endl;

		while (!inp.eof() && inp.peek() != '>') {
		    inp.get(buffer, 128, '>');
		}

		free(data);
		identifier = "Error reading";
		data = NULL;
		size = 0;

		delete e;
		return;
	    }
	}

	if (inp.eof() || inp.peek() == '>')
	    break;
    }

    data = (GeneticSymbols*)realloc(data, sizeof(GeneticSymbols) * size);
}


