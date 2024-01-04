# Vulkan FEM Solver

VulkanFEM is a FEM Solver implementation on the GPU using Vulkan. For the implementation of the GPU pipeline, [Kompute](https://github.com/KomputeProject/kompute/) is used. The problem created for solving is a basic Linear Elasticity Equation, where a voxel model is inputted to create FEM system and solve for displacement.

## Dependencies

- C++ compiler
  - On Windows: Visual Studio 2022 (for Visual C++ compiler)
  - On Linux: GCC
- CMake (tested with minimum 3.24)

## Submodules

The project depends on [Eigen](https://eigen.tuxfamily.org) and [Kompute](https://github.com/KomputeProject/kompute/) as submodules. Make sure to update the submodules after cloning the repository, with:

```bat
git submodule update --init
```

or download the projects manually into their folders.

## How to Run

This project is created to be easily used with Visual Studio Code with CMake, but a manual usage is also possible.

### With Visual Studio Code

Visual Studio Code has a very useful CMake extension, make sure to install CMake and enable the extension to use it for the project. With VS Code ready, just open the project folder with VS Code and it should automatically ask for building platform.

To build for the platform of choice, simply select build platform and configuration from the CMake bar on the bottom and hit build.

To run the project, Run and Debug window can be used. The current running profile must be matched manually to the CMake build configuration. The running profiles can be edited on [.vscode/launch.json](.vscode/launch.json).

### On Command Line/With Visual Studio

To build and run the project from the command line, CMake command line tools need to be used. In the project folder run

```bat
mkdir build
cd build
cmake -S .. --preset="Windows"
```

where preset is one of the configure presets according to the build platform, found in [CMakePresets.json](CMakePresets.json). After that

```bat
cd ..
cmake --build --preset "Windows Release"
```

with the corresponding build preset from [CMakePresets.json](CMakePresets.json). Alternatively on Windows, the solution created by CMake can be opened with Visual Studio.

To run the project,

```bat
cd .\build\Windows\src\Release\
.\VulkanFEM.exe data\armadillo_128.dat
```

or

```bash
cd build/Linux\ Release/src/
./VulkanFEM data/armadillo_128.dat
```

depending on the configurations selected.

### Known Build Issues

- A patch to Kompute is done as a step of CMake configuration, however this can fail in certain cases and cause the build to fail. In such cases, the patch might need to be applied manually. The patch can be found in [external/patches](external/patches/).
- In some cases during building with VS Code, compilation of the shaders can be done incompletely, as they consist of two steps. A second build can be done to ensure correct shaders are used if they are modified.

## Additional Information

- Data for the program can be found in [data](data) folder, where there is a sample data and additional zipped data.
- The output of the program is saved to the data folder found in the build folder, same folder as the executable.
- For different types of solver runs, definitions in [FEMDefines.hpp](include/FEMDefines.hpp) and [main.cpp](src/main.cpp) can be changed.
