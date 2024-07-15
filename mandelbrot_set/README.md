# Mandelbrot set in quad-precision
Implemented using [Normalized Iteration Count](https://en.m.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set#:~:text=Coloring,-algorithms)

**Pseudocode:**
```python
for each pixel (Px, Py) on the screen do
    x0:= scaled x coordinate of pixel (scaled to lie in the Mandelbrot X scale (-2.5, 1))
    y0:= scaled y coordinate of pixel (scaled to lie in the Mandelbrot Y scale (-1, 1))
    x:= 0.0
    y:= 0.0
    iteration:= 0
    max_iteration:= 1000
    // Here N = 2^8 is chosen as a reasonable bailout radius.

    while x*x + y*y ≤ (1 << 16) and iteration < max_iteration do
        xtemp:= x*x - y*y + x0
        y:= 2*x*y + y0
        x:= xtemp
        iteration:= iteration + 1

    // Used to avoid floating point issues with points inside the set.
    if iteration < max_iteration then
        // sqrt of inner term removed using log simplification rules.
        log_zn:= log(x*x + y*y) / 2
        nu:= log(log_zn / log(2)) / log(2)
        // Rearranging the potential function.
        // Dividing log_zn by log(2) instead of log(N = 1<<8)
        // because we want the entire palette to range from the
        // center to radius 2, NOT our bailout radius.
        iteration:= iteration + 1 - nu

    color1:= palette[floor(iteration)]
    color2:= palette[floor(iteration) + 1]
    // iteration % 1 = fractional part of iteration.
    color:= linear_interpolate(color1, color2, iteration % 1)
    plot(Px, Py, color)
```
## Usage
```
python setup.py build_ext --inplace
```
## Interesting areas
```
width = 800
height = 800
max_iter = 1000
center_r = -0.75
center_i = 0.0
zoom = 1.0
```
![image](https://github.com/user-attachments/assets/ae12fc00-6afd-4c27-b055-3983ae1ff30e)

```
width = 800
height = 800
max_iter = 1024
center_r = -0.7450
center_i = 0.10
zoom = 50
```
![image](https://github.com/user-attachments/assets/b8df9dbe-778b-4dd3-9591-bc2ab82032fd)
