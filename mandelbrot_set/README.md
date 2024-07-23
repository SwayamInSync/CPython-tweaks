# Comparing Mandelbrot set for `quad-precision`, `double` and `long double` data types
Theoretical reference is taken from [Plotting algorithms for the Mandelbrot set(Wikipedia)](https://en.m.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set)

| double | long double | __float128 |
|--------|-------------|------------|
| ![Mandelbrot set using double](path_to_double_image) | ![Mandelbrot set using long double](path_to_long_double_image) | ![Mandelbrot set using __float128](path_to_float128_image) |

## Installation
- Compile and build [Sleef](https://sleef.org/)

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
