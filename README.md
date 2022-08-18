# MARSIS
Data for scheduling of MARSIS Radar operations

&nbsp;

---

**Legend of data files:**

0.1. ORBIT INFO: main informations about spacecraft orbits from 04/01/2004 (orbit 1) to 01/01/2026 (orbit 27764)
0.2. RegionOfInterest: (x,y) coordinates defining the contour of the areas with higher observation preference
1. Ephemeris time: time in seconds of the data sampling. The delta time is of 1.866 seconds.
2. Local true solar time: 
3. Mars solar longitude: 
4. Mars sun distance: (remove ?)
5. Orbit number: id of the orbit travelled by the radar in the corresponding ephemeris time
6. Roughness: to upload
7. Solar zenith angle: angle between the sun's rays and the vertical direction in the corresponding ephemeris time. From this information it is possible to compute the sun elevation angle (= 90 - solar zenith angle)
8. Spacecraft altitude: altitude of the spacecraft with respect to Mars' ground in the corresponding ephemeris time
9. Sub_sc_latitude: latitude of the spacecraft in the corresponding ephemeris time
10. Sub_sc_longitude: longitude of the spacecraft in the corresponding ephemeris time
11. X: x-position of the radar in the corresponding ephemeris time
12. Y: y-position of the radar in the corresponding ephemeris time

