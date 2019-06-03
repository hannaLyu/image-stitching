#include "opencv2/opencv.hpp"
#include "opencv2/core/core.hpp"  
#include "opencv2/features2d/features2d.hpp"  
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"  
#include "opencv2/xfeatures2d.hpp"
#include <iostream>  
#include <vector> 

using namespace std;
using namespace cv;

int main() {
	

	/*translation image*/
	//Mat img_1 = imread("M:\\Documents\\KITTI\\000000.png");
	//Mat img_2 = imread("M:\\Documents\\KITTI\\000001.png");
	/*spin image*/
	/*Mat img_1 = imread("M:\\Documents\\KITTI\\bark.120.tiff");
	Mat img_2 = imread("M:\\Documents\\KITTI\\bark.000.tiff");*/
	/*zoom + spin image*/
	//Mat img_1 = imread("M:\\Documents\\KITTI\\Tiles_perspective_undistort.png");
	//Mat img_2 = imread("M:\\Documents\\KITTI\\Tiles_perspective_distort.png");
	/*image blur*/
	//Mat img_1 = imread("M:\\Documents\\KITTI\\norway.jpg");
	//Mat img_2 = imread("M:\\Documents\\KITTI\\norwayblur.jpg");
	/*light change*/
	Mat img_1 = imread("M:\\Documents\\KITTI\\plane1.jpg");
	Mat img_2 = imread("M:\\Documents\\KITTI\\planelighter.jpg");
	if (!img_1.data || !img_2.data)
	{
		cout << "error reading images " << endl;
		waitKey(0);
		return 0;
	}
    


	/*feature detect*/
	vector<KeyPoint> keyPoints_1, keyPoints_2;
	Mat descriptors_1, descriptors_2;
	//define feature point type
	int n = 3;



	/*-----------------FAST featrue Point----------------*/
	if (n == 1) {
		Ptr<FastFeatureDetector> detector_FAST = FastFeatureDetector::create();
		detector_FAST->detect(img_1, keyPoints_1);
		detector_FAST->detect(img_2, keyPoints_2);
		Ptr<ORB> descriptor_ORB = ORB::create();
		descriptor_ORB->compute(img_1, keyPoints_1, descriptors_1);
		descriptor_ORB->compute(img_2, keyPoints_2, descriptors_2);
		cout << "FAST feature point " << endl;
	}
	
	/*-----------------SURF featrue Point----------------*/
	else if(n==2){
		cv::Ptr<cv::xfeatures2d::SURF> surf = cv::xfeatures2d::SURF::create();
		surf->detect(img_1, keyPoints_1);
		surf->detect(img_2, keyPoints_2);
		surf->detectAndCompute(img_1, cv::Mat(), keyPoints_1, descriptors_1);
		surf->detectAndCompute(img_2, cv::Mat(), keyPoints_2, descriptors_2);
		cout << "SURF feature point " << endl;
	}
	/*-----------------SIFT featrue Point----------------*/
	else if(n==3){
		cv::Ptr<cv::xfeatures2d::SIFT> sift = cv::xfeatures2d::SIFT::create();
		sift->detect(img_1, keyPoints_1);
		sift->detect(img_2, keyPoints_2);
		sift->detectAndCompute(img_1, cv::Mat(), keyPoints_1, descriptors_1);
		sift->detectAndCompute(img_2, cv::Mat(), keyPoints_2, descriptors_2);
		cout << "SIFT feature point " << endl;
		
	}


	/*feature matching*/
	BFMatcher matcher;
	std::vector< DMatch > matches;
	matcher.match(descriptors_1, descriptors_2, matches);

	double max_dist = 0; double min_dist = 100;
	//Task 1:  Quick calculation of max and min distances between keypoints  
	for (int i = 0; i < descriptors_1.rows; i++)
	{
		double dist = matches[i].distance;
		if (dist < min_dist) min_dist = dist;
		if (dist > max_dist) max_dist = dist;
	}
	cout << "-- Max dist :" << max_dist << endl;
	cout << "-- Min dist :" << min_dist << endl;

	//-- Draw only "good" matches (i.e. whose distance is less than 0.6*max_dist )  
	//-- PS.- radiusMatch can also be used here.  
	vector< DMatch > good_matches;
	for (int i = 0; i < descriptors_1.rows; i++)
	{
		if (matches[i].distance < 0.6*max_dist)
		{
			good_matches.push_back(matches[i]);
		}
	}


	vector<KeyPoint> m_LeftKey;
	vector<KeyPoint> m_RightKey;
	vector<DMatch> m_Matches;


	int ptCount = (int)matches.size();
	Mat p1(ptCount, 2, CV_32F);
	Mat p2(ptCount, 2, CV_32F);

	// transfer Keypoint to Mat
	Point2f pt;
	for (int i = 0; i<ptCount; i++)
	{
		pt = keyPoints_1[matches[i].queryIdx].pt;
		p1.at<float>(i, 0) = pt.x;
		p1.at<float>(i, 1) = pt.y;

		pt = keyPoints_2[matches[i].trainIdx].pt;
		p2.at<float>(i, 0) = pt.x;
		p2.at<float>(i, 1) = pt.y;
	}


	// Task2: caculate F matrix based on RANSAC
	Mat m_Fundamental;
	vector<uchar> m_RANSACStatus;

	m_Fundamental = findFundamentalMat(p1, p2, m_RANSACStatus, FM_RANSAC);

	// outliers
	int OutlinerCount = 0;
	for (int i = 0; i<ptCount; i++)
	{
		if (m_RANSACStatus[i] == 0)
		{
			OutlinerCount++;
		}
	}

	// inliers
	vector<Point2f> m_LeftInlier;
	vector<Point2f> m_RightInlier;
	vector<DMatch> m_InlierMatches;
	int InlinerCount = ptCount - OutlinerCount;
	m_InlierMatches.resize(InlinerCount);
	m_LeftInlier.resize(InlinerCount);
	m_RightInlier.resize(InlinerCount);
	InlinerCount = 0;
	for (int i = 0; i<ptCount; i++)
	{
		if (m_RANSACStatus[i] != 0)
		{
			m_LeftInlier[InlinerCount].x = p1.at<float>(i, 0);
			m_LeftInlier[InlinerCount].y = p1.at<float>(i, 1);
			m_RightInlier[InlinerCount].x = p2.at<float>(i, 0);
			m_RightInlier[InlinerCount].y = p2.at<float>(i, 1);
			m_InlierMatches[InlinerCount].queryIdx = InlinerCount;
			m_InlierMatches[InlinerCount].trainIdx = InlinerCount;
			InlinerCount++;
		}
	}

	// 把内点转换为drawMatches可以使用的格式
	vector<KeyPoint> key1(InlinerCount);
	vector<KeyPoint> key2(InlinerCount);
	KeyPoint::convert(m_LeftInlier, key1);
	KeyPoint::convert(m_RightInlier, key2);

	// 显示计算F过后的内点匹配
	//Mat m_matLeftImage;
	//Mat m_matRightImage;
	// 以上两个变量保存的是左右两幅图像
	Mat OutImage;
	drawMatches(img_1, key1, img_2, key2, m_InlierMatches, OutImage);
	if (n == 1) {
		cout << "-- FAST InlinerCount :" << InlinerCount << endl;
	}
	else if (n == 2) {
		cout << "-- SURF InlinerCount :" << InlinerCount << endl;
	}
	else if (n == 3) {
		cout << "-- SIFT InlinerCount :" << InlinerCount << endl;
	}
	//stereoRectifyUncalibrated();

	Mat img_matches;
	drawMatches(img_1, keyPoints_1, img_2, keyPoints_2,
		good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
		vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);

	imwrite("FASTResult.jpg", img_matches);
	imshow("Match", img_matches);

	imwrite("FmatrixResult.jpg", OutImage);
	imshow("Match2", OutImage);
	waitKey(0);
	//system("pause");
	return 0;
}