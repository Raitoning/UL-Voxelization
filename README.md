# UL-Voxelization
Research initiation on voxelization techniques at the Université de Lorraine (France).

The goal of this project is to learn about different voxelization techniques, how they can be used and how to speed them up.

The project is a C++ implementation of different voxelization techniques, and uses [DGTal](https://github.com/DGtal-team/DGtal) as a graphic library to see and interact with the output.

# Compiling the project

To compile the project, you will first need to install [DGTal](https://github.com/DGtal-team/DGtal) on your machine. You can get it here: https://github.com/DGtal-team/DGtal.

## Compiling DGTal on Ubuntu 16.04
We ran into issues compiling DGTal on our machines running Ubuntu 16.04. Here is how we managed to get it installed (please note that you will need `cmake`, you can get it by typing `sudo apt install cmake` or something like this in a terminal. You might also need `boost`, `qt5` and other dependences. I can't help you much about it, it was a nightmare to compile.):

1. Download/Clone the latest version of DGTal on your computer.
2. Make a `build` folder in the cloned repository.
3. Run `cmake .. -DWITH_QGLVIEWER=true -DWITH_QT5=true` in the `build` directory.
4. Run `make -j [# of threads of your CPU]` (Will take some time, almost 5 minutes on an overclocked Intel Core i5 3570k @4.2Ghz x4!).
5. Run `sudo make install`.

At this point, DGTal should be installed as `libDGtal` on your computer (supposedly in `/usr/local/lib/libDGtal.so`, at least it's there on my desktop and laptop.) If you still have issues compiling DGTal, you might want to check [DGTal](https://github.com/DGtal-team/DGtal)'s repo for additionnal info/help. It took us 2 weeks to get it compiled. Have fun.

## Compiling the project

To compile the project, you will need to:

1. Download/Clone the latest version of the repo on your computer.
2. Make a `build` folder in the cloned repository.
3. Run `cmake ..` in the `build` directory.
4. Run `make`.
5. Run `./display3D file.off` to launch the executable.

If everything went fine, a window should open with a 3D bunny inside (or the file you chose). You can interact with it using the mouse and keyboard. Use the H key to open an help window with some keyboard/mouse actions.

In the near future, you will be able to choose what kind of voxelisation you want between Gaussian (1 voxel each unit) and a Resolution based one (128x128x128, 256x256x256, ...).

# Roadmap

The roadmap will be updated as the project progress and the tasks are assigned by our teacher.

- [X] Installing DGtal on our machines.
- [ ] Perform a Gaussian voxelization of a 3D model (1 voxel every unit.)
- [X] Perform a resolution based voxelization (ie 512x512x512 for example.)
- [ ] (Bonus) Implement the voxelization on a GPU (API to be chosen).
