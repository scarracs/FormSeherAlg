#ifndef M_PI
#define M_PI 3.141592654
#endif

#include "edl/edl.h"
#include "mathutil.h"

#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

#include <vector>
#include <iostream>


EDL::EDL(int sobelKernelSize, double sobelScale, double sobelDelta, int gaussianKernelSize,
         int minAnchorThreshold, int anchorStepping, int anchorThreshold, double angleTolerance, unsigned int minLineLength)
    : sobelKernelSize(sobelKernelSize),
      sobelScale(sobelScale),
      sobelDelta(sobelDelta),
      gaussianKernelSize(gaussianKernelSize),
      minAnchorThreshold(minAnchorThreshold),
      anchorThreshold(anchorThreshold),
      anchorStepping(anchorStepping),
      angleTolerance(angleTolerance),
      minLineLength(minLineLength)
{
}

EDL::~EDL(){}

std::vector<Line> EDL::calculate(cv::InputArray _image)
{
    image = _image.getMat();
    std::vector<Line> result;

    if(!image.data)
        return result;

    //Definitions for grad, angle and anchor calculation
    gradientMagnitudes = cv::Mat::zeros(image.rows, image.cols, CV_8U);
    dx = cv::Mat::zeros(image.rows, image.cols, CV_16S);
    dy = cv::Mat::zeros(image.rows, image.cols, CV_16S);
    std::vector<cv::Point> anchors;

    // ####
    // use a filter algorithm to suppress noise
    // ####

    cv::GaussianBlur(image, image, cv::Size(gaussianKernelSize, gaussianKernelSize), 0, 0);

    // ####
    // create gradient and direction output including the anchorpoints
    // ####

    calcGrad();
    findAnchors(anchors);

    // ####
    // run the routing algorithm
    // ####

    routeAnchors(angleTolerance, anchors, result);

    // Save result
    return result;
}

void EDL::findAnchors(std::vector<cv::Point> &anchors)
{
    bool direction;
    bool isAnchor = false;
    int nRows = gradientMagnitudes.rows-1;
    int nCols = gradientMagnitudes.cols-1;
    short center;

    for(int row = 1; row < nRows; ++row)
    {
        uchar* gradMag = gradientMagnitudes.ptr<uchar>(row);
        short* gradX = dx.ptr<short>(row);
        short* gradY = dy.ptr<short>(row);

        for(int column = 1; column < nCols; column += anchorStepping)
        {
            if ((center = gradMag[column]) > minAnchorThreshold)  // check if the current Point might be an anchor
            {
                cv::Vec2s vec(gradX[column] ,gradY[column]);
                direction = getOrientation(vec);
                if (direction == HORIZONTAL)
                {
                    if (center - (gradientMagnitudes.at<uchar>(row, column - 1)) >= anchorThreshold && center - (gradientMagnitudes.at<uchar>(row, column + 1) >= anchorThreshold ))
                    {
                        isAnchor = true;
                    }
                }
                if (direction == VERTICAL)
                {
                    if (center - (gradientMagnitudes.at<uchar>(row - 1, column)) >= anchorThreshold && center - (gradientMagnitudes.at<uchar>(row + 1, column) >= anchorThreshold ))
                    {
                        isAnchor = true;
                    }
                }
            }

            if(isAnchor) // save the Point
            {
                anchors.push_back(cv::Point(column, row));
                isAnchor = false;
            }
        }
    }
}

void EDL::calcGrad()
{
    int nRows = gradientMagnitudes.rows-1;
    int nCols = gradientMagnitudes.cols-1;

    short sobelX[3][3] =   {{-1,0,1},
                            {-2,0,2},
                            {-1,0,1}};

    short sobelY[3][3] =   {{-1,-2,-1},
                            {0,0,0},
                            {1,2,1}};

    short Gx; // sobel value in X-Direction
    short Gy; // sobel value in Y-Direction
    short gradientMagnitude; // length of the vector[Gx,Gy]

    for(int row = 1; row < nRows; ++row)
    {
        uchar* gradMag = gradientMagnitudes.ptr<uchar>(row);
        short* gradDx = dx.ptr<short>(row);
        short* gradDy = dy.ptr<short>(row);

        for(int column = 1; column < nCols; ++column)
        {
            Gx =    (sobelX[0][0] * image.at<uchar>(row-1, column-1)) + (sobelX[0][1] * image.at<uchar>(row, column-1)) + (sobelX[0][2] * image.at<uchar>(row+1,column-1)) +
                    (sobelX[1][0] * image.at<uchar>(row-1, column))   + (sobelX[1][1] * image.at<uchar>(row, column))   + (sobelX[1][2] * image.at<uchar>(row+1,column)) +
                    (sobelX[2][0] * image.at<uchar>(row-1, column+1)) + (sobelX[2][1] * image.at<uchar>(row, column+1)) + (sobelX[2][2] * image.at<uchar>(row+1,column+1));

            Gy =    (sobelY[0][0] * image.at<uchar>(row-1, column-1)) + (sobelY[0][1] * image.at<uchar>(row, column-1)) + (sobelY[0][2] * image.at<uchar>(row+1,column-1)) +
                    (sobelY[1][0] * image.at<uchar>(row-1, column))   + (sobelY[1][1] * image.at<uchar>(row, column))   + (sobelY[1][2] * image.at<uchar>(row+1,column)) +
                    (sobelY[2][0] * image.at<uchar>(row-1, column+1)) + (sobelY[2][1] * image.at<uchar>(row, column+1)) + (sobelY[2][2] * image.at<uchar>(row+1,column+1));

            gradientMagnitude = math::sqrtFast(Gx*Gx + Gy*Gy);
            gradMag[column] = gradientMagnitude;
            gradDx[column] = Gx;
            gradDy[column] = Gy;
        }
    }
}

void EDL::routeAnchors(double angleTolerance, std::vector<cv::Point>& anchorPoints, std::vector<Line> &result)
{
    cv::Mat_<uchar> edgels = cv::Mat::zeros(gradientMagnitudes.rows, gradientMagnitudes.cols, CV_8U);
    std::vector<std::list<cv::Point*>*> lineSegments;

    // Iterate all anchor points
    for(auto anchorPoint : anchorPoints)
    {
        // Is the pixel already part of an edge?
        if(edgels(anchorPoint) != 0)
        {
            continue;
        }

        walkFromAnchor(anchorPoint, lineSegments);
    }

    // Create result
    for(auto lineSegment : lineSegments)
    {
        if(lineSegment->size() >= minLineLength)
        {
            cv::Point* start = lineSegment->front();
            cv::Point* end = lineSegment->back();

            result.push_back(Line(*start, *end));
        }
    }

    // Free lineSegments
    for(auto lineSegment : lineSegments)
    {
        for(auto point : *lineSegment)
        {
            delete point;
        }
        delete lineSegment;
    }
}

void EDL::walkFromAnchor(cv::Point& anchorPoint, std::vector<std::list<cv::Point*>*>& lineSegments)
{
    //save the results to
    std::list<cv::Point*>* currentLineSegment = new std::list<cv::Point*>;
    double currentGradientAngle;
    double lineSegmentAngle;

    // two points to work with

    cv::Point *currentPoint = new cv::Point(anchorPoint);
    cv::Point *nextPoint = new cv::Point(anchorPoint);

    // Directions to walk to

    cv::Point *left = new cv::Point(0, -1);
    cv::Point *right = new cv::Point(0, 1);
    cv::Point *up = new cv::Point(1, 0);
    cv::Point *down = new cv::Point(-1, 0);

    int mainDirection; // VERTICAL or HORIZONTAL
    cv::Point *subDirection; // left,rigt or up,down

    // ####
    // set defaults before we start walking
    // ####

    bool stopWalk = false;
    cv::Vec2s vec(dx.at<short>(*currentPoint), dy.at<short>(*currentPoint));
    mainDirection = getOrientation(vec);

    if (mainDirection == HORIZONTAL)
    {
        subDirection = left;
    }

    if (mainDirection == VERTICAL)
    {
        subDirection = up;
    }

    // ####
    // Start the walking...
    // ####

    do
    {

        // ####
        // Check if a new segment begins
        // ####

        currentGradientAngle = fabs(std::atan2(dx.at<short>(*currentPoint), dy.at<short>(*currentPoint)));

        if (currentLineSegment->empty())
        {
            lineSegmentAngle = currentGradientAngle;
            currentLineSegment->push_front(currentPoint);
        }
        else
        {
            if(isAligned(lineSegmentAngle, currentGradientAngle, angleTolerance))
            {
                if(subDirection == left || subDirection == up)
                {
                    currentLineSegment->push_front(currentPoint);
                }
                else
                {
                    currentLineSegment->push_back(currentPoint);
                }
            }
            else
            {
                if(currentLineSegment->size() > 30) // if the line is to small forget about it
                {
                 lineSegments.push_back(currentLineSegment);
                }
                currentLineSegment = new std::list<cv::Point*>;
            }
        }

        // ####
        // Find next point
        // ####

        *nextPoint = getNextPoint(*currentPoint, mainDirection, *subDirection);

        // ####
        // Perform various checks
        // ####

        if(gradientMagnitudes.at<uchar>(*nextPoint) <= 0) // if its useless to keep walking into this direction
        {
            if (subDirection == right || subDirection == down) //if both directions already have been walked
            {
                stopWalk = true;
            }

            if (subDirection == left) // change from left to right
            {
                subDirection = right;
                *currentPoint = anchorPoint;
            }

            if (subDirection == up) // change from up to down
            {
                subDirection = down;
                *currentPoint = anchorPoint;
            }
        }

        else  // if we keep walking switch the points
        {
            currentPoint = nextPoint;
        }

        // will the next given point still be there?
        // TODO: Find better solution
        if( nextPoint->x+1 >= gradientMagnitudes.cols || nextPoint->y+1 >= gradientMagnitudes.rows || nextPoint->x < 1 || nextPoint->y < 1 )
        {
            stopWalk = true;
        }
    } while(!stopWalk);
}

bool EDL::getOrientation(cv::Vec2s &v1)
{
    bool orientation = VERTICAL;

    if(abs(v1[0]) >= abs(v1[1]))
    {
        orientation = HORIZONTAL;
    }
    return orientation;

}

cv::Point EDL::getNextPoint(cv::Point& currentPoint, int mainDirection, cv::Point& subDirection)
{
    int nRows;
    int nCols;
    int startRow = currentPoint.y;
    int startColumn = currentPoint.x;

    if (mainDirection == HORIZONTAL)
    {
        startRow -= 1;
        startColumn += subDirection.x;
        nRows = 3;
        nCols = 1;
    }
    else // VERTICAL
    {
        startRow += subDirection.y;
        startColumn -= 1;
        nRows = 1;
        nCols = 3;
    }

    cv::Point nextPoint;
    uchar currentMag = 0;
    uchar debugmag = 0 ;

    for(int i = 0 ; i < nRows; ++i)
    {
        uchar *gradMag = gradientMagnitudes.ptr<uchar>(startRow + i);

        for(int j = 0 ; j < nCols; ++j)
        {
            debugmag = gradMag[startColumn + j];
            if (currentMag < gradMag[startColumn + j])
            {
                currentMag = gradMag[startColumn + j];
                nextPoint = cv::Point(startColumn + j, startRow + i);
            }
        }
    }
    return nextPoint;
}

bool EDL::isAligned(double compare, double angle, double tolerance)
{
    if(compare >= M_PI)
        compare -= M_PI;
    if(angle >= M_PI)
        angle -= M_PI;

    return ((angle - compare) <= tolerance) && ((angle - compare) >= -tolerance);
}

bool EDL::isOutOfBounds(cv::Point *point)
{
    return (point->x < 0) || (point->x > gradientMagnitudes.cols)
            || (point->y < 0) || (point->y > gradientMagnitudes.rows);
}

double EDL::getAngleBetweenVectors(cv::Vec2s &v1, cv::Vec2s &v2)
{
    double cross = v1[0]*v2[1] - v1[1]*v2[0];
    double dot = v1.dot(v2);
    return atan2(cross, dot);
}
