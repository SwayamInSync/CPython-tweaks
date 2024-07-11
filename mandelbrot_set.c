/*
Z(t+1) = Z(t)^2 + C
Range in real line is [-2, 0.25]
*/


#include <stdio.h>
#include <math.h>

#define WIDTH 100
#define HEIGHT 80
#define THREHSOLD 4.0
#define MAX_ITER 1000

double get_abs_value(double r, double i)
{
    return sqrt(r*r + i*i);
}

int mandelbrot(double cr, double ci) 
{
    double z_real = 0.0, z_img = 0.0;
    double z_real2, z_img2;

    for (int n = 0; n < MAX_ITER; n++) {
        z_real2 = z_real * z_real;
        z_img2 = z_img * z_img;

        if (get_abs_value(z_real2, z_img2) > THREHSOLD)
            return n;

        z_img = 2.0 * z_real * z_img + ci;
        z_real = z_real2 - z_img2 + cr;
    }
    return MAX_ITER;
}

int main() {
    for (int y = 0; y < HEIGHT; y++) 
    {
        for (int x = 0; x < WIDTH; x++) {
            double cr = (x - WIDTH/2.0) * 4.0 / WIDTH; // linear transformation to keeping in the range
            double ci = (y - HEIGHT/2.0) * 4.0 / HEIGHT;
            
            int value = mandelbrot(cr, ci);

            if(value == MAX_ITER)
                printf("*");
            else
                printf(" ");
        }
        putchar('\n');
    }
    return 0;
}
