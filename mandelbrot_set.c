/*
Z(t+1) = Z(t)^2 + C
Range in real line is [-2, 0.25]
*/


#include <stdio.h>

#define WIDTH 80
#define HEIGHT 80
#define THREHSOLD 4.0
#define MAX_ITER 1000


int mandelbrot(double cr, double ci) {
    double zr = 0.0, zi = 0.0;
    double zr2, zi2;
    int n;
    for (n = 0; n < MAX_ITER; n++) {
        zr2 = zr * zr;
        zi2 = zi * zi;
        if (zr2 + zi2 > THREHSOLD)
            return n;
        zi = 2.0 * zr * zi + ci;
        zr = zr2 - zi2 + cr;
    }
    return MAX_ITER;
}

int main() {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            double cr = (x - WIDTH/2.0) * 4.0 / WIDTH; // linear transformation to keeping in the range
            double ci = (y - HEIGHT/2.0) * 4.0 / HEIGHT;
            
            int value = mandelbrot(cr, ci);
            
            // if (value == MAX_ITER)
            //     putchar(' ');
            // else if (value > MAX_ITER/2)
            //     putchar('.');
            // else if (value > MAX_ITER/4)
            //     putchar('*');
            // else if (value > MAX_ITER/8)
            //     putchar('+');
            // else
            //     putchar('#');

            if(value == MAX_ITER)
                printf("*");
            else
                printf(" ");
        }
        putchar('\n');
    }
    return 0;
}