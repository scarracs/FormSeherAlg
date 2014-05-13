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
        std::cout << test_gradientMagnitudes << std::endl << test_dx << std::endl << test_dy << std::endl;
     }

//    void findAnchorsTest()
//    {
//        //Create an image
//        cv::Mat_<uchar> image = (cv::Mat_<uchar>(4,5) <<    0,       0,       0,      0,       0,
//                                                            1,      11,     255,     21,       0,
//                                                            0,      12,     254,     22,       0,
//                                                            0,       0,       0,      0,       0);
//        //Check the image
//        QVERIFY(image(cv::Point(0, 1)) == 1);
//        QVERIFY(image(cv::Point(1, 1)) == 11);
//        QVERIFY(image(cv::Point(2, 1)) == 255);
//        QVERIFY(image(cv::Point(3, 1)) == 21);

//        //get the variables ready and set them
//        edl->image = image;
//        edl->gradientMagnitudes = cv::Mat::zeros(image.rows, image.cols, CV_8U);
//        edl->dx = cv::Mat::zeros(image.rows, image.cols, CV_16S);
//        edl->dy = cv::Mat::zeros(image.rows, image.cols, CV_16S);

//        //call the method
//        edl->calcGrad();

//        cv::Mat_<uchar> test_gradientMagnitudes = edl->gradientMagnitudes;
//        cv::Mat_<uchar> test_dx = edl->dx;
//        cv::Mat_<uchar> test_dy = edl->dy;
//        std::cout << test_gradientMagnitudes << std::endl << test_dx << std::endl << test_dy << std::endl;
//     }

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

        void getNextPointTest()
        {
            // HORIZONTAL
            cv::Point *left = new cv::Point(-1, 0);
            cv::Point *right = new cv::Point(1, 0);
            // VERTICAL
            cv::Point *up = new cv::Point(0, 1);
            cv::Point *down = new cv::Point(0, -1);
            // start/end
            cv::Point *currentPoint = new cv::Point(1, 1);
            cv::Point *nextPoint;
            // test 1
            cv::Mat_<uchar> gradientMagnitudes = (cv::Mat_<uchar>(3,3) <<        1        ,2      ,3
                                                                                ,11       ,255    ,13
                                                                                ,21       ,40     ,80); //walk over this
            edl->gradientMagnitudes = gradientMagnitudes;
            *nextPoint = edl->getNextPoint(*currentPoint, HORIZONTAL, *down);
            std::cout << *nextPoint << std::endl;
            QVERIFY(gradientMagnitudes(*nextPoint) == 80);
        }

//    void findNextPointTest()
//    {
//        cv::Mat_<uchar> magnitudes = cv::Mat::zeros(5, 5, CV_8U);
//        magnitudes(1, 1) = 150;
//        magnitudes(1, 2) = 140;
//        magnitudes(1, 3) = 143;
//        magnitudes(2, 2) = 155;
//        magnitudes(3, 2) = 140;
//        magnitudes(4,  3) = 120;

//        int mainDirection = HORIZONTAL;
//        int subDirection = -1;
//        cv::Point currentPoint(2, 1);

//        cv::Point* nextPoint = edl->findNextPoint(&currentPoint, mainDirection, subDirection, magnitudes);
//        QVERIFY(1 == nextPoint->x);
//        QVERIFY(1 == nextPoint->y);
//        QVERIFY(150 == magnitudes(*nextPoint));
//        delete nextPoint;

//        subDirection = +1;
//        nextPoint = edl->findNextPoint(&currentPoint, mainDirection, subDirection, magnitudes);
//        QVERIFY(3 == nextPoint->x);
//        QVERIFY(1 == nextPoint->y);
//        QVERIFY(143 == magnitudes(*nextPoint));
//        delete nextPoint;

//        mainDirection = VERTICAL;
//        subDirection = -1;
//        currentPoint.x = 2;
//        currentPoint.y = 3;

//        nextPoint = edl->findNextPoint(&currentPoint, mainDirection, subDirection, magnitudes);
//        QVERIFY(2 == nextPoint->x);
//        QVERIFY(2 == nextPoint->y);
//        QVERIFY(155 == magnitudes(*nextPoint));
//        delete nextPoint;

//        subDirection = +1;
//        nextPoint = edl->findNextPoint(&currentPoint, mainDirection, subDirection, magnitudes);
//        QVERIFY(3 == nextPoint->x);
//        QVERIFY(4 == nextPoint->y);
//        QVERIFY(120 == magnitudes(*nextPoint));
//        delete nextPoint;
//    }

    void isAlignedTest()
    {
        double angleTolerance = 21.5 * M_PI / 180.0d;

        // test aligned
        QVERIFY(true  == edl->isAligned( 0,   0, angleTolerance));
        QVERIFY(true  == edl->isAligned( 21.5 * M_PI / 180.0d, 0, angleTolerance));
        QVERIFY(true  == edl->isAligned(180 * M_PI / 180.0d,  0, angleTolerance));
        QVERIFY(true  == edl->isAligned(90 * M_PI / 180.0d,  270 * M_PI / 180.0d, angleTolerance));
        // Test not aligned
        QVERIFY(false == edl->isAligned( 22 * M_PI / 180.0d,   0, angleTolerance));
        QVERIFY(false == edl->isAligned( 45 * M_PI / 180.0d, 23 * M_PI / 180.0d, angleTolerance));
        QVERIFY(false == edl->isAligned(100 * M_PI / 180.0d, 200 * M_PI / 180.0d, angleTolerance));
    }

    void isOutOfBoundsTest()
    {
        cv::Mat mat = cv::Mat::zeros(2, 2, CV_8U);

        cv::Point point(-1, 0);
        QVERIFY(true == edl->isOutOfBounds(&point, mat));

        point = cv::Point(0, 0);
        QVERIFY(false == edl->isOutOfBounds(&point, mat));

        point = cv::Point(2, 2);
        QVERIFY(false == edl->isOutOfBounds(&point, mat));

        point = cv::Point(2, 3);
        QVERIFY(true == edl->isOutOfBounds(&point, mat));

        point = cv::Point(3, 2);
        QVERIFY(true == edl->isOutOfBounds(&point, mat));

        point= cv::Point(1, -1);
        QVERIFY(true == edl->isOutOfBounds(&point, mat));
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
};

#endif // EDLTEST_H
