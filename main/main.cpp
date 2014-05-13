#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <edl/edl.h>
#include <line.h>
#include <iostream>

using namespace std;
using namespace cv;

int main( int argc, char** argv )
{
    if(argc < 2)
    {
        std::cout << "Call with path to an image file!" << std::endl;
        return 0;
    }

    cv::Mat input = cv::imread(argv[1]);
    if(!input.data)
    {
        std::cerr << "Invalid input file. Must be an image!" << std::endl;
        std::cerr << "Could not load image '" << argv[1] << "'!" << std::endl;
        return -1;
    }
    Mat test;
 //   EDL edl = EDL();
  /// Create window
    cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );
    cv::imshow( "Display window", test );  //or grad or whatever

  waitKey(0);

  return 0;
  }
