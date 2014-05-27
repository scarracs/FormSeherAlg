#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <edl/edl.h>
#include <line.h>
#include <iostream>
#include <time.h>

using namespace std;
using namespace cv;

double get_time()
{
    struct timespec ts;
    if (clock_gettime (CLOCK_REALTIME, &ts) != 0)
    puts ("WARNING: Cannot read time using 'clock_gettime'!");
    return (double) ts.tv_sec + (double) ts.tv_nsec * 1e-9;
}

int main( int argc, char** argv )
{
    double startTime, endTime;
    int n;
    int count = 100;
    EDL *edl = new EDL();
    std::vector<Line> result;

   // Bild laden (wird nicht mitgemessen)

    if(argc < 2)
    {
        std::cout << "Call with path to an image file!" << std::endl;
        return 0;
    }

    cv::Mat input = cv::imread(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
    if(!input.data)
    {
        std::cerr << "Invalid input file. Must be an image!" << std::endl;
        std::cerr << "Could not load image '" << argv[1] << "'!" << std::endl;
        return -1;
    }

    // Zeimtessung beginnt

    startTime = get_time();
    for (n = 0; n < count ; n++)
    {
        edl->calculate(input);
    }
    endTime = get_time();
    printf ("\nElapsed time: %.2lf seconds\n\n", (endTime-startTime) / count);

    // Ausgabe des Ergebnisses (wird nicht mitgemessen)

    result = edl->calculate(input);

    cv::RNG rng(0xFFFFFFFF);
    cv::Mat resultImage = cv::Mat::zeros(input.rows, input.cols, CV_8UC3);
    for(auto line : result)
           {
               cv::Scalar color(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
               cv::line(resultImage, line.getStart(), line.getEnd(), color);
           }

  /// Create window
    cv::namedWindow( "Display window", cv::WINDOW_AUTOSIZE );
    cv::imshow( "Display window", resultImage );  //or grad or whatever

  waitKey(0);

  return 0;
  }
