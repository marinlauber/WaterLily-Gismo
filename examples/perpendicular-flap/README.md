### Perpendicular Flap

To start the Fluid (`WaterLily`) participant

```bash
module load precice/3.1.1
preCICEJulia --project=../ Fluid.jl precice-config.xml
```

To start the Solid (`G+Smo`) participant

```bash
module load precice/3.1.1
sudo make -C ~/Workspace/gismo/build/ test-2d-vertical-beam-implicit-solid
sudo ./test-2d-vertical-beam-implicit-solid -c precice-config.xml
```