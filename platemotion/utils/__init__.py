'''
plate subpackage

This subpackage defines some classes that facilitate the processing of plate motion.

Class structure:

    Plate

        - attributes:
            - vertices: vertices of a closed spherical polygon in form of [[lat_0,lon_0],...,[lat_n,lon_n]]
            - lats: latitudes of the spherical polygon in degrees
            - lons: longitudes of the spherical polygon in degrees
            - orientation: vertices arrangement; it can be counterclockwise or clockwise

            - methods:
            - contains_points: determine if a single point or multiple points are inside a spherical polygon.
            - area: calculate the area of a spherical polygon.
            - perimeter: calculate the perimeter of a spherical polygon.
            - centroid: identify the location of the centroid of a spherical polygon.
            - inertia: calculate the inertia tensor of a spherical polygon.
            
'''