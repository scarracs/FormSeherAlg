#ifndef EDLTEST_H
#define EDLTEST_H

#include <QObject>
#include <QtTest/QtTest>

#include <opencv2/core/core.hpp>

#include <iostream>

// brute force everything public (don't do this at home kids!)
#undef private
#undef protected
#define private public
#define protected public

#include "edl/edl.h"

// Restore visibilities
#undef private
#undef protected
#define protected protected
#define private private

class EDLTest : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase()
    {
        edl = new EDL();
    }

    void calcGradTest()
    {
        //Create an image
        cv::Mat_<uchar> image = (cv::Mat_<uchar>(4,5) <<    0,       0,       0,      0,       0,
                                                            1,      11,     255,     21,       0,
                                                            0,      12,     254,     22,       0,
                                                            0,       0,       0,      0,       0);
        //Check the image
        QVERIFY(image(cv::Point(0, 1)) == 1);
        QVERIFY(image(cv::Point(1, 1)) == 11);
        QVERIFY(image(cv::Point(2, 1)) == 255);
        QVERIFY(image(cv::Point(3, 1)) == 21);

        //get the variables ready and set them
        edl->image = image;
        edl->gradientMagnitudes = cv::Mat::zeros(image.rows, image.cols, CV_8U);
        edl->dx = cv::Mat::zeros(image.rows, image.cols, CV_16S);
        edl->dy = cv::Mat::zeros(image.rows, image.cols, CV_16S);

        //call the method
        edl->calcGrad();

        cv::Mat_<uchar> test_gradientMagnitudes = edl->gradientMagnitudes;
        cv::Mat_<uchar> test_dx = edl->dx;
        cv::Mat_<uchar> test_dy = edl->dy;
        std::cout << "calcGradTest:" << std::endl;
        std::cout << test_gradientMagnitudes << std::endl << test_dx << std::endl << test_dy << std::endl;
    }

    void findAnchorsTest()
    {
        //Create an image
        cv::Mat_<uchar> image = (cv::Mat_<uchar>(4,5) <<    0,       0,       0,      0,       0,
                                                            1,      11,     255,     21,       0,
                                                            0,      12,     254,     22,       0,
                                                            0,       0,       0,      0,       0);
        //Check the image
        QVERIFY(image(cv::Point(0, 1)) == 1);
        QVERIFY(image(cv::Point(1, 1)) == 11);
        QVERIFY(image(cv::Point(2, 1)) == 255);
        QVERIFY(image(cv::Point(3, 1)) == 21);

        //get the variables ready and set them
        edl->image = image;
        edl->gradientMagnitudes = cv::Mat::zeros(image.rows, image.cols, CV_8U);
        edl->dx = cv::Mat::zeros(image.rows, image.cols, CV_16S);
        edl->dy = cv::Mat::zeros(image.rows, image.cols, CV_16S);
        std::vector<cv::Point> anchors;

        //call the method
        edl->calcGrad();
        edl->findAnchors(anchors);

        cv::Mat_<uchar> test_gradientMagnitudes = edl->gradientMagnitudes;
        std::cout << "findAnchorsTest:" << std::endl;
        std::cout << test_gradientMagnitudes << std::endl;
        std::cout << "anchors:" << anchors << std::endl;
     }

    void getOrientationTest()
    {
        cv::Vec2s v1 = cv::Vec2s(2,4);
        QVERIFY(edl->getOrientation(v1) == VERTICAL);

        v1 = cv::Vec2s(4,2);
        QVERIFY(edl->getOrientation(v1) == HORIZONTAL);

        v1 = cv::Vec2s(3,7);
        QVERIFY(edl->getOrientation(v1) == VERTICAL);

        v1 = cv::Vec2s(7,3);
        QVERIFY(edl->getOrientation(v1) == HORIZONTAL);

        v1 = cv::Vec2s(5,5);
        QVERIFY(edl->getOrientation(v1) == HORIZONTAL);
    }

    void getNextPointTest()
    {
        // HORIZONTAL
        cv::Point left = cv::Point(-1, 0);
        cv::Point right = cv::Point(1, 0);
        // VERTICAL
        cv::Point up = cv::Point(0, -1);
        cv::Point down = cv::Point(0, 1);
        // start/end
        cv::Point currentPoint = cv::Point(1, 1);
        cv::Point *nextPoint;


        cv::Mat_<uchar> gradientMagnitudes = (cv::Mat_<uchar>(5,4) <<        1        ,2      ,3    ,4
                                                                            ,11       ,255    ,13   ,14
                                                                            ,21       ,22     ,80   ,15
                                                                            ,31       ,32     ,33   ,16
                                                                            ,41       ,42     ,43   ,17
                                                                            ,51       ,52     ,53   ,18); //walk over this
        edl->gradientMagnitudes = gradientMagnitudes;

        //simple walk

        nextPoint = edl->getNextPoint(currentPoint, right);
        QVERIFY(gradientMagnitudes(*nextPoint) == 80);

        nextPoint = edl->getNextPoint(currentPoint, left);
        QVERIFY(gradientMagnitudes(*nextPoint) == 21);

        nextPoint = edl->getNextPoint(currentPoint, up);
        QVERIFY(gradientMagnitudes(*nextPoint) == 3);

        nextPoint = edl->getNextPoint(currentPoint, down);
        QVERIFY(gradientMagnitudes(*nextPoint) == 80);

        //corner walk!

        currentPoint = cv::Point(0, 3);

        nextPoint = edl->getNextPoint(currentPoint, right);
        QVERIFY(gradientMagnitudes(*nextPoint) == 42);

        nextPoint = edl->getNextPoint(currentPoint, left);
        QVERIFY(gradientMagnitudes(*nextPoint) == 31);

        nextPoint = edl->getNextPoint(currentPoint, up);
        QVERIFY(gradientMagnitudes(*nextPoint) == 22);

        nextPoint = edl->getNextPoint(currentPoint, down);
        QVERIFY(gradientMagnitudes(*nextPoint) == 42);
    }

    void walkFromAnchorTest()
    {
        //Create an image
        cv::Mat_<uchar> image = (cv::Mat_<uchar>(4,5) <<    0,       0,       0,      0,       0,
                                                            1,      11,     255,     21,       0,
                                                            0,      12,     254,     22,       0,
                                                            0,       0,       0,      0,       0);

        //get the variables ready and set them
        edl->image = image;
        edl->gradientMagnitudes = cv::Mat::zeros(image.rows, image.cols, CV_8U);
        edl->dx = cv::Mat::zeros(image.rows, image.cols, CV_16S);
        edl->dy = cv::Mat::zeros(image.rows, image.cols, CV_16S);
        std::vector<cv::Point> anchors;

        //call the previous methods
        edl->calcGrad();
        edl->findAnchors(anchors);

        //call walkFromAnchor
        std::vector<std::list<cv::Point*>*> lineSegments;
        edl->walkFromAnchor(anchors[0], lineSegments);

//        //show the result
//        std::cout << "walkFromAnchorTest:" << std::endl;
//        std::cout << "Line" << lineSegments[0] << std::endl;
    }

    void getAngleBetweenVectorsTest()
    {
        cv::Vec2s v1 = cv::Vec2s(2,4);
        cv::Vec2s v2 = cv::Vec2s(2,5);
        double angle1 = calcAngle(v1,v2);
        double angle2 = fabs(edl->getAngleBetweenVectors(v1, v2));
        QVERIFY((angle2 - angle1) > -0.01 && (angle2 - angle1) <= 0.01 );

        v1 = cv::Vec2s(3,7);
        v2 = cv::Vec2s(7,3);
        angle1 = calcAngle(v1,v2);
        angle2 = fabs(edl->getAngleBetweenVectors(v1, v2));
        QVERIFY((angle2 - angle1) > -0.01 && (angle2 - angle1) <= 0.01 );

        v1 = cv::Vec2s(1,4);
        v2 = cv::Vec2s(1,7);
        angle1 = calcAngle(v1,v2);
        angle2 = fabs(edl->getAngleBetweenVectors(v1, v2));
        QVERIFY((angle2 - angle1) > -0.01 && (angle2 - angle1) <= 0.01 );
    }

    void isOutOfBoundsTest()
    {
        edl->gradientMagnitudes = cv::Mat::zeros(2,3,CV_8U);

        cv::Point point(-1, 0);
        QVERIFY(true == edl->isOutOfBounds(point));
        QVERIFY(true == edl->isOutOfBounds(-1, 0));

        point = cv::Point(0, 0);
        QVERIFY(false == edl->isOutOfBounds(point));
        QVERIFY(false == edl->isOutOfBounds(0, 0));

        point = cv::Point(2, 2);
        QVERIFY(false == edl->isOutOfBounds(point));
        QVERIFY(false == edl->isOutOfBounds(2,2 ));

        point = cv::Point(2, 3);
        QVERIFY(true == edl->isOutOfBounds(point));
        QVERIFY(true == edl->isOutOfBounds(2, 3));

        point = cv::Point(3, 2);
        QVERIFY(false == edl->isOutOfBounds(point));
        QVERIFY(false == edl->isOutOfBounds(3, 2));

        point= cv::Point(1, -1);
        QVERIFY(true == edl->isOutOfBounds(point));
        QVERIFY(true == edl->isOutOfBounds(1, -1));
    }

    void cleanupTestCase()
    {
        delete edl;
    }

private:
    EDL* edl;

    double calcAngle(cv::Vec2s &v1, cv::Vec2s &v2)
    {
        return std::acos(v1.dot(v2) / ((std::sqrt(v1[0]*v1[0]+v1[1]*v1[1])) * (std::sqrt(v2[0]*v2[0]+v2[1]*v2[1]))));
    }

    double degreeToRadiant(double angle)
    {
        return angle * M_PI / 180.0d;
    }

    double radiantToDegree(double radiant)
    {
        return radiant * 180.0d / M_PI;
    }
};

#endif // EDLTEST_H
