#include<stdio.h>
#include<sleef.h>
#include<sleefquad.h>

#define WIDTH 100
#define HEIGHT 80
#define THRESHOLD 4.0Q
#define MAX_ITER 1000

Sleef_quad get_abs_value(Sleef_quad *r, Sleef_quad *i)
{
    Sleef_quad rr = Sleef_mulq1_u05(*r, *r);
    Sleef_quad ii = Sleef_mulq1_u05(*i, *i);
    return Sleef_sqrtq1_u05(Sleef_addq1_u05(rr, ii));
}

int mandelbrot(Sleef_quad *cr, Sleef_quad *ci)
{
    Sleef_quad z_real = 0.0Q;
    Sleef_quad z_img = 0.0Q;
    Sleef_quad z_real2, z_img2;

    for(int i = 0; i < MAX_ITER; i++)
    {
        z_real2 = Sleef_mulq1_u05(z_real, z_real);
        z_img2 = Sleef_mulq1_u05(z_img, z_img);

        if(get_abs_value(&z_real2, &z_img2) > THRESHOLD)
            return i;
        
        z_img = Sleef_addq1_u05(Sleef_mulq1_u05(2.0Q, Sleef_mulq1_u05(z_real, z_img)), *ci);
        z_real = Sleef_addq1_u05(Sleef_subq1_u05(z_real2, z_img2), *cr);
    }

    return MAX_ITER;
}

int main()
{
    for(int y = 0; y < HEIGHT; y++)
    {
        for(int x = 0; x < WIDTH; x++)
        {
            Sleef_quad cr = Sleef_mulq1_u05(Sleef_subq1_u05(x, WIDTH/2.0Q), 4.0Q/WIDTH);
            Sleef_quad ci = Sleef_mulq1_u05(Sleef_subq1_u05(y, HEIGHT/2.0Q), 4.0Q/HEIGHT);  

            int value = mandelbrot(&cr, &ci);
            if(value == MAX_ITER)
                printf("*");
            else
                printf(" ");
        }
        printf("\n");
    }

    return 0;
}
