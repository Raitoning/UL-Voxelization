#include <string>
#include <DGtal/base/Common.h>
#include <DGtal/io/Display3D.h>
#include <DGtal/io/readers/MeshReader.h>
#include <DGtal/io/viewers/Viewer3D.h>

using namespace DGtal;

// Function to split a string into multiple strings.
// Code from Cplusplus.com
// http://www.cplusplus.com/articles/2wA0RXSz/
const std::vector<std::string> Explode(const std::string &s, const char &c)
{
    std::string buff{""};
    std::vector<std::string> v;

    for (auto n : s)
    {
        if (n != c)
            buff += n;
        else if (n == c && buff != "")
        {
            v.push_back(buff);
            buff = "";
        }
    }
    if (buff != "")
        v.push_back(buff);

    return v;
}

int main(int argc, char **argv)
{
    std::string fileName;
    std::string line;
    std::vector<std::string> values;

    // Checking the arguments.
    // If there isn't a file name as argument, stop the execution.
    if (argc < 2)
    {
        DGtal::trace.info() << "Usage: Viewer fileName" << std::endl
                            << "This program can open voxels file (.vox) or an object file (.off)" << std::endl
                            << "Now exiting..." << std::endl;

        return 0;
    }

    std::string importFileName = argv[1];

    // qT Application hosting the viewer.
    QApplication application(argc, argv);

    // 3D Viewer.
    Viewer3D<> viewer;

    if (importFileName.find(".vox") != std::string::npos)
    {
        ifstream file(argv[1]);

        double x;
        double y;
        double z;

        std::vector<Z3i::RealPoint> voxels;

        if (file.is_open())
        {
            int i = 0;
            while (std::getline(file, line))
            {
                values = Explode(line, ' ');

                x = atof(values[0].c_str());
                y = atof(values[1].c_str());
                z = atof(values[2].c_str());

                voxels.emplace_back(Z3i::RealPoint(x, y, z));
            }

            file.close();

            std::string importFileName = argv[1];
            int dirSize = 0;
            for (int i = 0; i < importFileName.size(); i++)
            {
                if (importFileName.at(i) == '/')
                {
                    dirSize = i;
                }
            }

            importFileName = importFileName.substr(dirSize);

            std::cout << "Succesfully imported voxels from " << importFileName << std::endl;
        }
        else
        {
            std::cout << "Failed to import voxels." << std::endl;
        }

        for (uint i = 0; i < voxels.size(); i++)
        {
            viewer.addCube(voxels[i]);
        }
    }
    else if (importFileName.find(".off") != std::string::npos)
    {
        Mesh<Viewer3D<>::RealPoint> mesh;

        mesh << importFileName;

        viewer << mesh;
    }
    else
    {
        std::cout << "The file isn't a voxel (.vox) file or a mesh (.off) file." << std::endl
                  << "Now exiting..." << std::endl;

        return 0;
    }

    viewer << Viewer3D<>::updateDisplay;
    viewer.show();

    // Return the qT application.
    return application.exec();
}
