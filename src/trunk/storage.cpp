#include "storage.h"

template<class T>
void parse(T** data, size_t* size, istream& inp) {
    inp >> *size;
    *data = (T*)malloc(sizeof(T) * *size);

    size_t aSize = 0;
    for (size_t i=0; i<(*size); ++i) {
	char c;
	inp.get(c);
	if (c == ' ' || c == '\t' || c == '\n' || c == '\r')
	    continue;
	(*data)[aSize++] = fromChar(c);
    }

    *size = aSize;
}

template void parse<GeneticSymbols>(GeneticSymbols** data, size_t* size, istream& inp);

