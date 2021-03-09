# Welcome to the SphericalPolygon package

[![PyPI version shields.io](https://img.shields.io/pypi/v/sphericalpolygon.svg)](https://pypi.python.org/pypi/sphericalpolygon/) [![PyPI pyversions](https://img.shields.io/pypi/pyversions/sphericalpolygon.svg)](https://pypi.python.org/pypi/sphericalpolygon/) [![PyPI status](https://img.shields.io/pypi/status/sphericalpolygon.svg)](https://pypi.python.org/pypi/sphericalpolygon/) [![GitHub contributors](https://img.shields.io/github/contributors/lcx366/SphericalPolygon.svg)](https://GitHub.com/lcx366/SphericalPolygon/graphs/contributors/) [![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/lcx366/SphericalPolygon/graphs/commit-activity) [![GitHub license](https://img.shields.io/github/license/lcx366/SphericalPolygon.svg)](https://github.com/lcx366/SphericalPolygon/blob/master/LICENSE) [![Documentation Status](https://readthedocs.org/projects/pystmos/badge/?version=latest)](http://sphericalpolygon.readthedocs.io/?badge=latest)

The SphericalPolygon package is an archive of scientific routines for handling spherical polygons. Currently, operations on spherical polygons include:

1. calculate the area or mass(if the area density is given) 
2. calculate the perimeter
3. identify the centroid 
4. compute the geometrical or physical moment of inertia tensor
5. determine whether one or multiple points are inside the spherical polygon

## How to Install

On Linux, macOS and Windows architectures, the binary wheels can be installed using **pip** by executing one of the following commands:

```
pip install sphericalpolygon
pip install sphericalpolygon --upgrade # to upgrade a pre-existing installation
```

## How to use

### Create a spherical polygon

Spherical polygons can be created from a 2d array in form of `[[lat_0,lon_0],..,[lat_n,lon_n]]` with unit of degrees, or from a boundary file, such as those in [Plate boundaries for NNR-MORVEL56 model](http://geoscience.wisc.edu/~chuck/MORVEL/PltBoundaries.html). The spherical polygon accepts a latitude range of [-90,90] and a longitude range of [-180,180] or [0,360].


```python
from sphericalpolygon import Sphericalpolygon
# build a spherical polygon for Antarctica Plate
polygon = Sphericalpolygon.from_file('NnrMRVL_PltBndsLatLon/an',skiprows=1) 
print(polygon.orientation)
```

    Counterclockwise


### Calculate the area

Calculate the area(or the solid angle) of a spherical polygon over a unit sphere.


```python
print(polygon.area())
```

    1.4326235943514618


Calculate the area of the spherical polygon over the Earth with an averaged radius of 6371km.


```python
print(polygon.area(6371), ' km2')
```

    58149677.38285546  km2


Calculate the mass of the spherical polygon shell with a thickness of 100km and density of 3.1g/cm3 over the Earth with an averaged radius of 6371km.


```python
print(polygon.area(6371,100*3.1), ' Gt')
```

    18026399988.685192  Gt


### Calculate the perimeter

Calculate the perimeter of a spherical polygon over a unit sphere.


```python
print(polygon.perimeter())
```

    6.322665894174733


Calculate the perimeter of a spherical polygon over the Earth with an averaged radius of 6371km.


```python
print(polygon.perimeter(6371), ' km')
```

    40281.70441178723  km


### Identify the centroid

Identify the centroid of a spherical polygon over a unit sphere.


```python
print(polygon.centroid())
```

    [-83.61081032380656, 57.80052886741483, 0.13827778179537997]


Identify the centroid of a spherical polygon over the Earth with an averaged radius of 6371km.


```python
print(polygon.centroid(6371),' deg deg km')
```

    [-83.61081032380656, 57.80052886741483, 880.9677478183658]  deg deg km


It shows that the latitude of the centroid is close to the South Pole, and the centroid is located about 881km underground.

### Compute the moment of inertia tensor

Compute the geometrical moment of inertia tensor of a spherical polygon over a unit sphere. The tensor is symmetrical and has six independent components. The first three components are located diagonally, corresponding to $Q_{11}$, $Q_{22}$, and $Q_{33}$; the last three components correspond to $Q_{12}$, $Q_{13}$, and $Q_{23}$.


```python
print(polygon.inertia())
```

    [ 1.32669154  1.17471081  0.36384484 -0.05095381  0.05246122  0.08126929]


Compute the physical moment of inertia tensor of the spherical polygon shell with a thickness of 100km and density of 3.1g/cm3 over the Earth with an averaged radius of 6371km.


```python
print(polygon.inertia(6371,100*3.1)/1e12, ' Gt·Gm2')
```

    [677582.33535848 599961.08075046 185826.79201142 -26023.68204226
      26793.56591716  41506.73569238]  Gt·Gm2


### Points are inside a polygon？

Determine if a single point or multiple points are inside a given spherical polygon.

#### single point


```python
print(polygon.contains_points([75,152]))
```

    False


#### multiple points


```python
print(polygon.contains_points([[-85,130],[35,70]]))
```

    [True, False]


### Change log

- **1.2.1 — Feb 23,  2021**
  - Replace the function *create_polygon* for building a spherical polygon object from a 2d array with methods `from_array` and `from_file`.
- **1.2.0 — Mar 20,  2020**
  - Add the `perimeter()` method that may calculate the perimeter of a spherical polygon.
  - Add the `centroid()` method that may determaine the centroid location for a spherical polygon.

## Reference

Chunxiao, Li. "Inertia Tensor for MORVEL Tectonic Plates." *ASTRONOMICAL RESEARCH AND TECHNOLOGY* 13.1 (2016).
