#define PY_SSIZE_T_CLEAN

#include<Python.h>

#include <stdio.h>
#include <sleef.h>
#include <sleefquad.h>
#include <omp.h>
#include <math.h>

typedef struct {
    unsigned char r, g, b;
} Color;

typedef struct
{
    float smooth_iter;
    Sleef_quad final_z_real;
    Sleef_quad final_z_img;
} MandelbrotResult;

typedef struct
{
    Sleef_quad real;
    Sleef_quad img;
} ComplexQuad;

Sleef_quad get_abs_value(Sleef_quad *r, Sleef_quad *i) 
{
    Sleef_quad rr = Sleef_mulq1_u05(*r, *r);
    Sleef_quad ii = Sleef_mulq1_u05(*i, *i);
    return Sleef_addq1_u05(rr, ii);
}

Color get_color(double t, Sleef_quad *final_z_real, Sleef_quad *final_z_img) 
{
    Color color;
    t = fmod(t * 20, 1.0); // Adjust color cycle frequency for more bands
    const double epsilon = 1e-10;

    // Interior coloring
    if (fabs(t - 1.0) < epsilon) {
        Sleef_quad abs_z = Sleef_sqrtq1_u05(get_abs_value(final_z_real, final_z_img));
        double interior_t = Sleef_cast_to_doubleq1(Sleef_fmodq1(abs_z, Sleef_cast_from_doubleq1(1.0)));
        color.r = (unsigned char)(255 * interior_t);
        color.g = (unsigned char)(128 * (1 - interior_t));
        color.b = (unsigned char)(64 * interior_t);
        return color;
    }

    if (t < 0.16) 
    {
        // Dark blue to light blue
        double v = t / 0.16;
        color.r = 0;
        color.g = (unsigned char)(150 * v);
        color.b = (unsigned char)(255 * (0.5 + 0.5 * v));
    } 
    else if (t < 0.42) 
    {
        // Light blue to white
        double v = (t - 0.16) / 0.26;
        color.r = (unsigned char)(255 * v);
        color.g = (unsigned char)(150 + 105 * v);
        color.b = 255;
    } 
    else if (t < 0.6425) 
    {
        // White to yellow
        double v = (t - 0.42) / 0.2225;
        color.r = 255;
        color.g = 255;
        color.b = (unsigned char)(255 * (1 - v));
    } 
    else if (t < 0.8575) 
    {
        // Yellow to orange
        double v = (t - 0.6425) / 0.215;
        color.r = 255;
        color.g = (unsigned char)(255 * (1 - 0.4 * v));
        color.b = 0;
    } 
    else 
    {
        // Orange to dark blue (cycle completion)
        double v = (t - 0.8575) / 0.1425;
        color.r = (unsigned char)(255 * (1 - v));
        color.g = (unsigned char)(153 * (1 - v));
        color.b = (unsigned char)(128 * v);
    }

    return color;
}

MandelbrotResult mandelbrot(Sleef_quad *cr, Sleef_quad *ci, int MAX_ITER, Sleef_quad *RADIUS2) 
{
    Sleef_quad z_real = Sleef_cast_from_doubleq1(0.0);
    Sleef_quad z_img = Sleef_cast_from_doubleq1(0.0);
    Sleef_quad two = Sleef_cast_from_doubleq1(2.0);
    Sleef_quad z_real2, z_img2;

    int i;
    MandelbrotResult result;

    for (i = 0; i < MAX_ITER; i++) 
    {
        z_real2 = Sleef_addq1_u05(Sleef_subq1_u05(Sleef_mulq1_u05(z_real, z_real), Sleef_mulq1_u05(z_img, z_img)), *cr);
        z_img2 = Sleef_addq1_u05(Sleef_mulq1_u05(two, Sleef_mulq1_u05(z_real, z_img)), *ci);

        if (Sleef_icmpgtq1(get_abs_value(&z_real2, &z_img2), *RADIUS2)) 
        {
            Sleef_quad log_zn =  Sleef_divq1_u05(Sleef_logq1_u10(get_abs_value(&z_real2, &z_img2)), two);
            Sleef_quad nu =  Sleef_divq1_u05(Sleef_logq1_u10(Sleef_divq1_u05(log_zn, SLEEF_M_LN2q)), SLEEF_M_LN2q);
            
            result.smooth_iter = i + 1 - Sleef_cast_to_doubleq1(nu);
            result.final_z_real = z_real2;
            result.final_z_img = z_img2;
            return result;
        }
        
        z_real = z_real2;
        z_img = z_img2;
    }

    result.smooth_iter = MAX_ITER;
    result.final_z_real = z_real;
    result.final_z_img = z_img;
    return result;
}

ComplexQuad complex_quad_add(ComplexQuad a, ComplexQuad b)
{
    ComplexQuad res;
    res.real = Sleef_addq1_u05(a.real, b.real);
    res.img = Sleef_addq1_u05(a.img, b.img);
    return res;
}

ComplexQuad complex_quad_mul(ComplexQuad a, ComplexQuad b)
{
    ComplexQuad res;
    res.real = Sleef_subq1_u05(Sleef_mulq1_u05(a.real, b.real), Sleef_mulq1_u05(a.img, b.img));
    res.img = Sleef_addq1_u05(Sleef_mulq1_u05(a.real, b.img), Sleef_mulq1_u05(a.img, b.real));
    return res;
}

ComplexQuad complex_quad_pow(ComplexQuad z, int n)
{
    ComplexQuad res = {Sleef_cast_from_doubleq1(1.0), Sleef_cast_from_doubleq1(0.0)};
    for(int i=0; i < n; i++)
    {
        res = complex_quad_mul(res, z);
    }

    return res;
}

ComplexQuad complex_quad_derivative(ComplexQuad z, ComplexQuad c, int n)
{
    if(n == 0)
        return (ComplexQuad){Sleef_cast_from_doubleq1(1.0Q), Sleef_cast_from_doubleq1(0.0)};
    
    ComplexQuad prev = complex_quad_derivative(z, c, n-1);
    return complex_quad_mul(prev, (ComplexQuad){Sleef_cast_from_doubleq1(2.0), Sleef_cast_from_doubleq1(0.0)});
}

Sleef_quad estimate_distance(Sleef_quad cr, Sleef_quad ci, int period) 
{
    ComplexQuad c = {cr, ci};
    ComplexQuad z = c;
    for (int i = 0; i < period; i++) 
    {
        z = complex_quad_add(complex_quad_pow(z, 2), c);
    }
    
    ComplexQuad dz = complex_quad_derivative(z, c, period);
    ComplexQuad dc = {Sleef_cast_from_doubleq1(1.0), Sleef_cast_from_doubleq1(0.0)};
    for (int i = 0; i < period; i++) 
    {
        dc = complex_quad_add(complex_quad_mul(dc, (ComplexQuad){Sleef_cast_from_doubleq1(2.0), Sleef_cast_from_doubleq1(0.0)}), (ComplexQuad){Sleef_cast_from_doubleq1(1.0), Sleef_cast_from_doubleq1(0.0)});
    }
    
    Sleef_quad abs_dz = Sleef_sqrtq1_u05(get_abs_value(&dz.real, &dz.img));
    Sleef_quad abs_dc = Sleef_sqrtq1_u05(get_abs_value(&dc.real, &dc.img));
    
    return Sleef_divq1_u05(Sleef_mulq1_u05(abs_dc, Sleef_logq1_u10(abs_dz)), abs_dz);
}

void find_interesting_point(Sleef_quad *center_r, Sleef_quad *center_i, Sleef_quad zoom) 
{
    Sleef_quad step = Sleef_divq1_u05(Sleef_cast_from_doubleq1(1.0), zoom);
    Sleef_quad min_distance = Sleef_cast_from_doubleq1(INFINITY);
    Sleef_quad best_r = *center_r;
    Sleef_quad best_i = *center_i;
    
    for (int i = -10; i <= 10; i++) 
    {
        for (int j = -10; j <= 10; j++) 
        {
            Sleef_quad r = Sleef_addq1_u05(*center_r, Sleef_mulq1_u05(Sleef_cast_from_int64q1(i), step));
            Sleef_quad i = Sleef_addq1_u05(*center_i, Sleef_mulq1_u05(Sleef_cast_from_int64q1(j), step));
            Sleef_quad distance = estimate_distance(r, i, 100);  // Assuming a max period of 100
            if (Sleef_icmpltq1(distance, min_distance)) 
            {
                min_distance = distance;
                best_r = r;
                best_i = i;
            }
        }
    }
    
    *center_r = best_r;
    *center_i = best_i;
}

static PyObject* mandelbrot_set(PyObject* self, PyObject* args) 
{
    int width, height, max_iter;
    const char *center_r_str, *center_i_str, *zoom_str;

    if (!PyArg_ParseTuple(args, "iiisss", &width, &height, &max_iter, &center_r_str, &center_i_str, &zoom_str))
        return NULL;
    unsigned char* img = (unsigned char *)malloc(width * height * 3);
    
    if (!img) {
        PyErr_NoMemory();
        return NULL;
    }

    Sleef_quad RADIUS = Sleef_cast_from_doubleq1(2.0);
    Sleef_quad RADIUS2 = Sleef_mulq1_u05(RADIUS, RADIUS);
    Sleef_quad zoom_q = Sleef_strtoq(zoom_str, NULL);
    zoom_q = Sleef_divq1_u05(1.0Q, zoom_q);

    Sleef_quad center_r = Sleef_strtoq(center_r_str, NULL);
    Sleef_quad center_i = Sleef_strtoq(center_i_str, NULL);

    #pragma omp parallel for
    for (int y = 0; y < height; y++) 
    {
        for (int x = 0; x < width; x++) 
        {
            Sleef_quad cr = Sleef_mulq1_u05(Sleef_divq1_u05(Sleef_subq1_u05(Sleef_cast_from_int64q1(x), Sleef_cast_from_int64q1(width / 2)), Sleef_cast_from_int64q1(width / 2)), RADIUS);
            Sleef_quad ci = Sleef_mulq1_u05(Sleef_divq1_u05(Sleef_subq1_u05(Sleef_cast_from_int64q1(y), Sleef_cast_from_int64q1(height / 2)), Sleef_cast_from_int64q1(height / 2)), RADIUS);
            
            cr = Sleef_addq1_u05(Sleef_mulq1_u05(cr, zoom_q), center_r);
            ci = Sleef_addq1_u05(Sleef_mulq1_u05(ci, zoom_q), center_i);

            MandelbrotResult result = mandelbrot(&cr, &ci, max_iter, &RADIUS2);
            
            Color color;
            if (result.smooth_iter == max_iter) 
            {
                // Interior point
                Sleef_quad abs_z = Sleef_sqrtq1_u05(get_abs_value(&result.final_z_real, &result.final_z_img));
                double interior_t = Sleef_cast_to_doubleq1(Sleef_fmodq1(abs_z, Sleef_cast_from_doubleq1(1.0)));
                color = get_color(1.0, &result.final_z_real, &result.final_z_img);
            } 
            else 
            {
                double t = result.smooth_iter / (double)max_iter;
                color = get_color(t, &result.final_z_real, &result.final_z_img);
            }
            
            img[(y * width + x) * 3 + 0] = color.r;
            img[(y * width + x) * 3 + 1] = color.g;
            img[(y * width + x) * 3 + 2] = color.b;
        }
    }

    PyObject* result = PyBytes_FromStringAndSize((char*)img, width * height * 3);
    free(img);
    return result;
}

static PyMethodDef MandelbrotMethods[] = 
{
    {"mandelbrot_set", mandelbrot_set, METH_VARARGS, "Calculate Mandelbrot set"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef mandelbrotmodule = 
{
    PyModuleDef_HEAD_INIT,
    "mandelbrot",
    NULL,
    -1,
    MandelbrotMethods
};

PyMODINIT_FUNC PyInit_mandelbrot(void) 
{
    return PyModule_Create(&mandelbrotmodule);
}
