#ifndef EDL_H
#define EDL_H

#ifndef M_PI
#define M_PI 3.141592654
#endif

#include "algorithm.h"

#include "line.h"

#include "mathutil.h"

#include <opencv2/core/core.hpp>
#include <vector>
#include <list>

#define HORIZONTAL 0
#define VERTICAL   1

class EDL : public Algorithm
{
public:
    EDL(int gaussianKernelSize = 3, int minAnchorThreshold = 40, int anchorStepping = 1, int anchorThreshold = 30, double angleTolerance = 20 *  M_PI / 180.0, unsigned int minLineLength = 30);
    ~EDL();

    std::vector<Line> calculate(cv::InputArray _image);

private:
    void calcGrad();

    void findAnchors(std::vector<cv::Point> &anchors);

    void routeAnchors(std::vector<cv::Point>& anchorPoints, std::vector<Line> &result);

    void walkFromAnchor(cv::Point& anchorPoint, std::vector<std::list<cv::Point*>*>& lineSegments);

    bool getOrientation(cv::Vec2s &v1);

    cv::Point* getNextPoint(cv::Point& currentPoint, cv::Point& subDirection);

    bool isOutOfBounds(cv::Point &point);

    bool isOutOfBounds(int x, int y);

    cv::Vec2s* getSobelVector(cv::Point &point);

    cv::Vec2s* getSobelVector(int x, int y);

    double getAngleBetweenVectors(cv::Vec2s &v1, cv::Vec2s &v2);

    int gaussianKernelSize;
    int anchorThreshold;
    double angleTolerance;
    unsigned int minLineLength;
    int minAnchorThreshold;
    int anchorStepping;
    cv::Mat image;
    cv::Mat gradientMagnitudes;
    cv::Mat dx;
    cv::Mat dy;
};

#endif // EDL_H
