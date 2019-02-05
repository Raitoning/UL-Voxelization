#include "DGtal/io/readers/MeshReader.h"
#include "DGtal/io/Display3D.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/base/Common.h"
#include "DGtal/io/Color.h"

using namespace std;
using namespace DGtal;

int main( int argc, char** argv )
{
  QApplication application(argc,argv);
  Viewer3D<> viewer;
  viewer.show();
  std::string inputFilename = "../Mesh/bunny.off";
  if( argc > 1)
    inputFilename = argv[1];

  // Since the input points are not necessary integers we use the PointD3D from Display3D.
  Mesh<Viewer3D<>::RealPoint> anImportedMesh;
  anImportedMesh << inputFilename;
  trace.info()<< "importating done..."<< endl;
  viewer.setLineColor(DGtal::Color(150,0,0,254));
  viewer << anImportedMesh;
  viewer << Viewer3D<>::updateDisplay;
  return application.exec();
}