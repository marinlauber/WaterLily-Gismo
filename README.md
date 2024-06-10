# WaterLily-Gismo

WaterLily &amp; G+Smo fluid-structure interaction through the `PreCICE` coupling library

### Prerequisites

The installation will assume that you have a working `PreCICE` library (`v3.1.1`). If not, please follow the instructions on the [PreCICE website](https://precice.org/). You will also need the `Julia` bindings for `PreCICE`. You can find the instructions to install the bindings [here](https://github.com/precice/PreCICE.jl).

You should also install `G+Smo` on your machine. You can find the instructions on the [G+Smo website](https://gismo.github.io/).


### Installation

To run G+Smo and WaterLily coupled simulations, you will need to clone the `PreCICE` brach of the `WaterLily` repository. You can do this by running the following command:
```bash
git clone git@github.com:marinlauber/WaterLily.git
cd WaterLily
git checkout PreCICE
```
You will also need the `ParametricBodies.jl` package in the `PlanarBodies` branch. You can clone the repository and checkout the branch by running the following commands:
```bash
git clone https://github.com/marinlauber/ParametricBodies.jl
cd ParametricBodies.jl
git checkout PlanarBodies
```

### Running the simulation

To run the example in this repo, you will first need to "dev" `WaterLily` and `ParameticBodies` into the environment `examples/Projects.toml`
```bash
cd WaterLily-Gismo/examples/perpendicular-flap
preCICEjulia --project=../
] dev /path/to/ParametricBodies.jl
] dev /path/to/WaterLily
] instantiate
```
Where in my case the command `preCICEjulia` is an alias to the julia executable with the necessary libraries preloaded. You can create the alias by adding the following line to your `.bashrc` or `.bash_aliases` file:
```bash
alias preCICEJulia='LD_PRELOAD="/usr/lib/x86_64-linux-gnu/libstdc++.so.6:/usr/lib/x86_64-linux-gnu/libcurl.so.4" julia'
```
You will then be able to run the simulation by running the following command in to separate terminals:
```bash
preCICEjulia --project=../ Fluid.jl precice-config.xml
sudo /path/to/gismo/bin/Solid -c precice-config.xml
```
