#include "1Dproblem.h"

int main()
{
    omp_set_num_threads(10);
    SodTubeProblem();
    return 0;
}