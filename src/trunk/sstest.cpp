#include "util.h"
#include "storage.h"

int main(void) {
    srandom(time(NULL));
    size_t arr[3000];

    for (int i=0; i<3000; i++)
	if ((i % 10) == 0)
	    arr[i] = 1600;
	else
	    arr[i] = 40;

    size_t res[3000];
    
    for (size_t i=0; i<300; i++)
	res[stochasticSelect(arr, 3000)]++;

    for (int i=0; i<3000; i++) {
	if (res[i] == 0) {
	    if ((i%10) == 0)
		printf("%d:     *\n", i);
	    else
		printf("%d: \n", i);
	} else {
	    if ((i%10) == 0)
		printf("%d: %4d*\n", i, res[i]);
	    else
		printf("%d: %4d\n", i, res[i]);    
	}
    }

    //printf("%d %d %d %d %d\n", res[0], res[1], res[2], res[3], res[4]);
}

/*int main(void) {
    GISeq 
    }*/
