# hyperdraw
Drawing octagonal tilings of hyperbolic space

------

This is a Mathematica script that will generate and draw an octagonal tiling of the Poincare disk. The original intent was to get a natural pictoral representation of the Calyey graph of the fundamental group of a genus 2 Riemann surface. The corresponding Riemann surface is obtained by gluing together the geodesics with the same label/color in the disk.

The tiling is generated by first determining the base octagon whose center is at the center of the disk. The coordinates of its vertices determine the type of tiling. There are infinitely many different tilings that can occur, and they are parameterized by the interior angle of the largest octagon, which must be of the form pi/(2n) for n>=1. Included in the repository are some examples for pi/2 and pi/4. Once the base octagon is generated, it is hyperbolically reflected about each of its sides recursively. 

Included is also a python script (which is untested) that will generate the tiling in a faster manner. I won't provide details on it because it is not finished. 
