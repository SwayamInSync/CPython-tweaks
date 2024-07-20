#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include <stdio.h>
#include <sleef.h>
#include <sleefquad.h>
#include <omp.h>
#include <math.h>

typedef struct {
    unsigned char r, g, b;
} Color;

typedef struct {
    float smooth_iter;
    Sleef_quad final_z_real;
    Sleef_quad final_z_img;
    Sleef_quad min_dist;
} MandelbrotResult;

typedef struct {
    Sleef_quad real;
    Sleef_quad img;
} ComplexQuad;

Sleef_quad get_abs_value(Sleef_quad *r, Sleef_quad *i) 
{
    Sleef_quad rr = Sleef_mulq1_u05(*r, *r);
    Sleef_quad ii = Sleef_mulq1_u05(*i, *i);
    return Sleef_addq1_u05(rr, ii);
}

Color get_color(double t, double interior_t, double orbit_trap) 
{
    Color color;
    const double epsilon = 1e-10;

    // Smooth transition between interior and exterior
    double blend = fmax(0, fmin(1, (t - 0.99) * 100));
    t = t * (1 - blend) + interior_t * blend;

    // Apply orbit trap influence
    t = t * 0.7 + orbit_trap * 0.3;

    // HSV to RGB conversion with smooth transitions
    double h = fmod(t * 10, 6);
    double s = 0.8;
    double v = 1 - pow(1 - t, 4);

    double c = v * s;
    double x = c * (1 - fabs(fmod(h, 2) - 1));
    double m = v - c;

    double r, g, b;
    if (h < 1) { r = c; g = x; b = 0; }
    else if (h < 2) { r = x; g = c; b = 0; }
    else if (h < 3) { r = 0; g = c; b = x; }
    else if (h < 4) { r = 0; g = x; b = c; }
    else if (h < 5) { r = x; g = 0; b = c; }
    else { r = c; g = 0; b = x; }

    color.r = (unsigned char)((r + m) * 255);
    color.g = (unsigned char)((g + m) * 255);
    color.b = (unsigned char)((b + m) * 255);

    return color;
}

MandelbrotResult mandelbrot(Sleef_quad *cr, Sleef_quad *ci, int MAX_ITER, Sleef_quad *RADIUS2) 
{
    Sleef_quad z_real = Sleef_cast_from_doubleq1(0.0);
    Sleef_quad z_img = Sleef_cast_from_doubleq1(0.0);
    Sleef_quad two = Sleef_cast_from_doubleq1(2.0);
    Sleef_quad z_real2, z_img2;
    Sleef_quad min_dist = Sleef_cast_from_doubleq1(INFINITY);

    int i;
    MandelbrotResult result;

    for (i = 0; i < MAX_ITER; i++) 
    {
        z_real2 = Sleef_addq1_u05(Sleef_subq1_u05(Sleef_mulq1_u05(z_real, z_real), Sleef_mulq1_u05(z_img, z_img)), *cr);
        z_img2 = Sleef_addq1_u05(Sleef_mulq1_u05(two, Sleef_mulq1_u05(z_real, z_img)), *ci);

        // Orbit trapping
        Sleef_quad dist = Sleef_sqrtq1_u05(get_abs_value(&z_real2, &z_img2));
        if (Sleef_icmpltq1(dist, min_dist)) {
            min_dist = dist;
        }

        if (Sleef_icmpgtq1(get_abs_value(&z_real2, &z_img2), *RADIUS2)) 
        {
            Sleef_quad log_zn = Sleef_logq1_u10(get_abs_value(&z_real2, &z_img2));
            Sleef_quad nu = Sleef_logq1_u10(log_zn / SLEEF_M_LN2q) / SLEEF_M_LN2q;
            
            result.smooth_iter = i + 1 - Sleef_cast_to_doubleq1(nu);
            result.final_z_real = z_real2;
            result.final_z_img = z_img2;
            result.min_dist = min_dist;
            return result;
        }
        
        z_real = z_real2;
        z_img = z_img2;
    }

    result.smooth_iter = MAX_ITER;
    result.final_z_real = z_real;
    result.final_z_img = z_img;
    result.min_dist = min_dist;
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

typedef struct {
    ComplexQuad z;
    ComplexQuad dz;
    ComplexQuad dc;
    ComplexQuad dzdz;
} Derivatives;

Derivatives iterate_and_compute_derivatives(ComplexQuad c, int max_iter) {
    ComplexQuad z = {Sleef_cast_from_doubleq1(0.0), Sleef_cast_from_doubleq1(0.0)};
    ComplexQuad dz = {Sleef_cast_from_doubleq1(1.0), Sleef_cast_from_doubleq1(0.0)};
    ComplexQuad dc = {Sleef_cast_from_doubleq1(0.0), Sleef_cast_from_doubleq1(0.0)};
    ComplexQuad dzdz = {Sleef_cast_from_doubleq1(0.0), Sleef_cast_from_doubleq1(0.0)};
    ComplexQuad two = {Sleef_cast_from_doubleq1(2.0), Sleef_cast_from_doubleq1(0.0)};

    for (int i = 0; i < max_iter; i++) {
        dzdz = complex_quad_add(complex_quad_mul(two, complex_quad_mul(z, dzdz)), complex_quad_mul(dz, dz));
        dz = complex_quad_add(complex_quad_mul(two, complex_quad_mul(z, dz)), dc);
        z = complex_quad_add(complex_quad_mul(z, z), c);
        dc.real = Sleef_cast_from_doubleq1(1.0);
        dc.img = Sleef_cast_from_doubleq1(0.0);
    }

    Derivatives result = {z, dz, dc, dzdz};
    return result;
}

Sleef_quad estimate_interior_distance(ComplexQuad c, int max_iter) {
    Derivatives d = iterate_and_compute_derivatives(c, max_iter);
    
    Sleef_quad dz_abs_sq = get_abs_value(&d.dz.real, &d.dz.img);
    Sleef_quad one = Sleef_cast_from_doubleq1(1.0);
    Sleef_quad numerator = Sleef_subq1_u05(one, dz_abs_sq);
    
    ComplexQuad denominator_term1 = complex_quad_mul(d.dc, d.dz);
    ComplexQuad denominator_term2 = complex_quad_mul(d.dzdz, complex_quad_mul(d.z, d.dc));
    ComplexQuad denominator = complex_quad_add(denominator_term1, denominator_term2);
    
    Sleef_quad denominator_abs = Sleef_sqrtq1_u05(get_abs_value(&denominator.real, &denominator.img));
    
    return Sleef_divq1_u05(numerator, denominator_abs);
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

    Sleef_quad RADIUS = Sleef_cast_from_doubleq1(2.0);  // Increased from 2.0
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
            
         
            double t, interior_t, orbit_trap;
            if (result.smooth_iter >= max_iter) 
            {
                ComplexQuad c = {cr, ci};
                Sleef_quad distance = estimate_interior_distance(c, max_iter);
                interior_t = Sleef_cast_to_doubleq1(Sleef_fmodq1(distance, Sleef_cast_from_doubleq1(1.0)));
                t = 1.0;
                orbit_trap = 0.0;
            } 
            else 
            {
                Sleef_quad distance = estimate_distance(cr, ci, 100);
                t = result.smooth_iter / (double)max_iter;
                t = 1 - pow(1 - t, 4);  // Non-linear smoothing
                orbit_trap = 1.0 - Sleef_cast_to_doubleq1(Sleef_divq1_u05(result.min_dist, RADIUS));
                interior_t = 0.0;
            }
            
            Color color = get_color(t, interior_t, orbit_trap);
            
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
