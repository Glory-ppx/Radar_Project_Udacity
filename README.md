# SFND - Udacity - Radar Target Generation and Detection

Here's my solution to the "Radar Project" from the Udacity's Sensor Fusion Nanodegree.
In this project, we simulate a FMCW radar with the following specifications:
- Frequency of operation: 77 GHz

```matlab
Frequency_of_operation = 77e9;
```

- Max Range: 200m

```matlab
MaxRange = 200;
```

- Range Resolution: 1m
```matlab
RangeResolution = 1;
```

- Max Velocity: 100m/s

```matlab
MaxVelocity = 100;
```

- Velocity Resolution: 1m/s

```matlab
VelResolution = 1;
```

- Speed of light: 300 000 000 m/s

```matlab
c = 3e8;
```

 - Bandwidth
 ```matlab
Bandwidth = c/(2*VelResolution);
 ```

 - Chirp time
 ```matlab
 Tchirp = 5.5*2*MaxRange/c;
 ```

As an example we set the target's initial position and velocity as 100m and 50 m/s, respectively. More over, we assume that its velocity remains contant.
```matlab
target_initial_position = 100;
target_velocity = 50;
```
