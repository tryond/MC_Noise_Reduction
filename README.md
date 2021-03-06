[![Generic badge](https://img.shields.io/badge/build-passing-<COLOR>.svg)](https://shields.io/)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

# Monte Carlo Rendering Noise Reduction

MC Rendering Noise Reduction is a program which implements several different filtering techniques for use with reducing noise generated by Monte Carlo Path Tracing renderings of 3D scenes.

## Description

Monte Carlo Path Tracing casts a ray from the viewer perspective out through each pixel in the frame. Each ray reflects off the objects in the scene a specified amount of times and returns a value for its respective pixel. As one can imagine, this operation can take a very long time to complete. This is compounded by the fact that Monte Carlo Renderings generally require many samples-per-pixel to produce clean results. 

The purpose of this project is to attain the results of a scene generated using high samples-per-pixel by applying post-processing filters to their low sample counterparts. 

### Goals:
1. Reduce overall noise
2. Maintain or improve upon time complexity
3. Maintain or improve upon visual results
4. Gain insights that could help us move towards real-time path tracing

## Motivation

MC Rendering Noise Reduction was created as my senior project in **CS114: Projects in Advanced 3D Computer Graphics** at the University of California, Irvine.

## Instructions

### Compile

Console compile command:

```console
mc-noise-reduce <username>$ g++ -std=c++11 -O3 -fopenmp -o simplept RNG.cpp Vec.cpp SimpleFilter.cpp ObjectFilter.cpp BrdfFilter.cpp NormalFilter.cpp DofFilter.cpp simplept.cpp
```

### Run

Console run command:

```console
mc-noise-reduce <username>$ simplept X a b c d e f g h i j k l m n
```

**Replace each letter in the run command with an intenger >= 0**

* **X**: Samples-per-pixel

**The following lower-case letters refer to filter types. The number designated refers to the size of the filter. Setting a filter to zero will turn the filter off entirely.**

* **a**) Simple mean parameter

* **b**) Simple median parameter

* **c**) Simple gaussian parameter

* **d**) Object mean parameter

* **e**) Object median parameter

* **f**) Object gaussian parameter

* **g**) Brdf mean parameter

* **h**) Brdf gaussian parameter

* **i**) Normal mean parameter

* **j**) Normal gaussian parameter

* **k**) Dof mean level parameter

* **l**) Dof mean focus level parameter

* **m**) Dof gray level parameter

* **n**) Dof gray focus level parameter

**A full write-up describing each filter can be found [here](https://tryond.github.io/noise_red.html).**

### Example

The following command will generate a **64** sample-per-pixel (**X**) rendering which applies both a **7** pixel wide simple median filter (**b**) and a **3** pixel wide simple gaussian filter (**c**). 

```console
mc-noise-reduce <username>$ simplept 64 0 7 3 0 0 0 0 0 0 0 0 0 0 0
```

![Example Rendering Image](res/combo_example.jpg?raw=true "Image that shows example rendering")

## Results

![Combo Filter gif](res/combo_16_anim.gif?raw=true "Animation that shows a scene rendered with a combo filter at 16spp")

* 16 samples-per-pixel
* 7 pixel wide simple median filter and a 3 pixel wide simple gaussian filter
* Render Time: 412.1 seconds

![No Filter gif](res/none_32_anim.gif?raw=true "Animation that shows a scene rendered with no filter at 32spp")

* 32 samples-per-pixel
* No filter
* Render Time: 815.7 seconds

After experimenting with the different types of filters and fine tuning their parameters, I found that the combination of the simple median and the simple gaussian worked very well. I was able to maintain the visual quality of the rendering with the higher sampling rate while cutting the time taken in half.

Please feel free to experiment and find out what combinations work well for you!

## Author

**Danny Tryon** - [tryond](https://github.com/tryond?tab=repositories)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

