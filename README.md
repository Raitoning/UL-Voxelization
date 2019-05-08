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

# Running the executables

If the compilation was succesful, you can run multiple programs.

## Voxelizer
If you want to voxelize a mesh, run the `Voxelizer` program with the file name and some of the following argments:

    --resolution="x y z"        Define the resolution used to voxelize

    --normalize=size            Define the minimal size of the mesh

    --threaded=# of threads     Define the number of threads to use. Default is the number of threads of the CPU.

    --origin="x y z"            Define the origin of a non orthogonal plane.

    --direction="x y z"         Define the direction of a non orthogonal plane.

    --delta=value               Define the delta of a non orthogonal plane.

    --export                    Define if the result must be exported to a .vox file.

Example:

    ./Voxelizer ../Mesh/bunny.off --normalize=8 --export --threaded=4

    ./Voxelizer ../Mesh/bunny.off --normalize=8 --resolution="8 8 8"

## Viewer

If you want to visualize the content of a mesh file (.off) or the result of a voxelization via a voxels file (.vox), run the `Viewer` program with the name of the file.

Examples:

    ./Viewer ../Mesh/bunny.off
    ./Viewer bunny_gaussian_10_10_8.vox

# Other infos

If everything went fine, a window should open with a 3D bunny inside (or the file you chose). You can interact with it using the mouse and keyboard. Use the H key to open an help window with some keyboard/mouse actions.

You can choose what kind of voxelisation you want between Gaussian (1 voxel each unit) and a Resolution based one (128x128x128, 256x256x256, ...).

As it is a CPU based voxelization, it is generally a bad idea to go over a normalization or a resolution of 64.
The time needed for voxelizing a mesh will mainly depends on the number of face of the mesh, as a ray must be tested for intersection for every faces.
For example, you can safely voxelize the `kitten_poisson2` or `blobby` meshes up to a definition of 32, or even 64, as the meshes contains a few faces. But for the Stanford bunny mesh (`bunny`), it is generally not advised to go over a definition of 16, as the mesh contains over 60k (60 000) faces.
Finding a decent definition for a mesh without having to wait too long is trial and error. I recommend you starting at a low definition, such as 4, then multiplying each time per 2 (4, 8, 16, ...) until you find a sweetspot between definition and computation time.

By default, the `Voxelizer` software will use as much threads (including SMT technologies such as Intel Hyper-Threading) as there are on the CPU of the machine, ex: 1 on a mono core CPU, 2 on a dual core CPU, 4 on a quad core CPU, ... If you want more control on the number of threads used, you can use the `--threaded` argment to specify how much threads to use. You can completely disable multithreading by specifying `1` as the number of threads (`--threaded=1`).

# Known bugs and limitations

* At the moment, only the Gaussian mode (without the `--resolution` argment) will compute voxels. As a result, the exported voxels file from a resolution based voxelization with be empty.

* If the voxelization process takes too long with the `Voxelizer` software, the data won't be passed to `DGtal` and `Qt` for display and result in a crash. However, the crash should happen after the export process, so you should be able to export the result to a voxel file via the `--export` option, and use the `Viewer` software to display the result.
**Note:** if you use multiple threads, you can perform more complex voxelization in the same time as a simple voxelization on a single thread, allowing you to display more data in `DGtal` without crashing.

* You can use the `--origin`, `--direction` and `--delta` argments of the `Voxelizer` software, but they will be ignored as there is an ongoing issue with the generation of a custom to send rays. Either a resolution based or a Gaussian voxelization will be performed, based on the argments you passed to the software.

# Roadmap

The roadmap will be updated as the project progress and the tasks are assigned by our teacher.

- [X] Installing DGtal on our machines.
- [X] Perform a Gaussian voxelization of a 3D model (1 voxel every unit.)
- [ ] Perform a resolution based voxelization (ie 512x512x512 for example.)
- [ ] Perfom a voxelization from a custom plan (ongoing.)
- [X] Export voxelization results to a .vox file.
- [X] Display voxels from a .vox file.
- [X] (Bonus) Implement the voxelization with multithreading support (OpenMP).
