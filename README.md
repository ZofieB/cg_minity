# minity - A minimal grpahics engine


## Getting Started

minity was developed under Windows using Microsoft Visual Studio 2019. It uses [CMake](https://cmake.org/) as its build system and does not contain anything else platform-specific, so it should in principle work using other compilers and operating systems, but this has not been tested. Should there be any issues, we are grateful for any input that helps us make the software work on as many platforms as possible.

## Prerequisites

The project uses [CMake](https://cmake.org/) and relies on the following libraries: 

- [GLFW](https://www.glfw.org/) 3.3 or higher (https://github.com/glfw/glfw.git) for windowing and input support
- [glm](https://glm.g-truc.net/) 0.9.9.5 or higher (https://github.com/g-truc/glm.git) for its math funtionality
- [glbinding](https://github.com/cginternals/glbinding) 3.1.0 or higher (https://github.com/cginternals/glbinding.git) for OpenGL API binding
- [globjects](https://github.com/cginternals/globjects) 2.0.0 or higher (https://github.com/cginternals/globjects.git) for additional OpenGL wrapping
- [Dear ImGui](https://github.com/ocornut/imgui) 1.71 or higher (https://github.com/ocornut/imgui.git) for GUI elements
- [tinyfiledialogs](https://sourceforge.net/projects/tinyfiledialogs/) 3.3.9 or higher (https://git.code.sf.net/p/tinyfiledialogs/code) for dialog functionality
- [stb](https://github.com/nothings/stb/) (stb_image.h 2.22 or higher and stb_image_write.h 1.13 or higher) or higher (https://github.com/nothings/stb.git) for image loading and saving

The release package includes all these libraries in compiled form (where applicable), so if you are using it you can skip the remainder of the section. If you want or need to compile the dependencies yourself, on Windows we provide a set of batch scripts to simplify the process of retrieving and building them. These scripts require that the Git and CMake executables properly set up and in your system path.

- ```./fetch-libs.cmd``` retrieves the libraries from their respective git repositories.  
- ```./build-libs.cmd``` builds them and installs them into the ```./lib``` folder.  
- ```./copy-libs.cmd``` copies DLLs required for execution into the ```./bin/Debug``` and ```./bin/Release``` folders. 

Under Windows, after opening a command prompt and running these three scripts from the project root folder, all dependencies should be available in the ```./lib``` folder.

```
./fetch-libs.cmd
./build-libs.cmd
./copy-libs.cmd
```

On other platforms, as Dear ImGui, and tinyfiledialogs, and stb do not use CMake, place them as subfolders of the ```./lib``` folder (using the folder names ```imgui```, ```tinyfd```, and ```stb```):

```
cd lib
git clone https://github.com/ocornut/imgui.git  
git clone http://git.code.sf.net/p/tinyfiledialogs/code tinyfd
git clone https://github.com/nothings/stb.git
cd ..
```

The other libraries use CMake, so just follow their respective instructions and make sure that your CMake installation can find them.

## Building

If you are using Visual Studio 2019, you can use its integrated CMake support to build and run the project. When running the program from within the IDE, make sure the current working directory is set to the project root in ```./.vs/launch.vs.json``` file by including ```"currentDir": "${projectDir}“```, so that all resources can be found by the executable.

When instead building from the command line, run the following commands from the project root folder:

```
mkdir build
cd build
cmake ..
```

After this, you can compile the debug or release versions of the project using 

```
cmake --build --config Debug
```

and/or

```
cmake --build --config Release
```

After building, the executables will be available in the ```./bin/Debug``` and ```./bin/Release``` folders.

## Running

As mentioned above, the program requires that the current working directory is set to the project root folder. On Windows, the program further assumes that the DLLs for glfw, glbinding, and globjects are located in the same folder as the executable (or, alternatively, are in the system path). If you followed the steps above, this should already be the case. Otherwise, you can run the copy-libs.cmd script to copy them over from the ```./lib``` folder. When running the program from outside the IDE, you can use the ```./minity-debug.cmd``` and ```./minity-release.cmd``` scripts in the project root folder to run the software.

## Usage

After starting the program, a file dialog will pop up and ask you for a Wavefront OBJ File file. Some basic usage instructions are displayed in the console window.

