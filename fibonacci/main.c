#include<stdio.h>
#include<sleef.h>
#include<sleefquad.h>

int main()
{
    Sleef_quad x = Sleef_cast_from_doubleq1(0.1);
    Sleef_quad result = Sleef_sinq1_u10(x);

    Sleef_printf("Result %Qe\n", result);

    return 0;
}

/*
cmd: export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib && gcc main.c -o main -lsleef -lsleefquad && ./main
*/