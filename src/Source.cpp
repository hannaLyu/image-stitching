#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include "opencv2/core.hpp"
#include "opencv2/features2d.hpp"
#include "opencv2/highgui.hpp"
using namespace cv;
//using namespace std;


// Global variable declaration
const Mat img_KITTI_0 = imread("KITTI\000000.png", IMREAD_GRAYSCALE);
const Mat img_KITTI_1 = imread("KITTI\000001.png", IMREAD_GRAYSCALE);
		
// Function 'mean_time_detectors'
template<class type> float mean_time_detectors(Ptr<type> &detector) {
	std::vector<KeyPoint> keypoints_KITTI;
	auto start = std::chrono::high_resolution_clock::now();
	detector->detect(img_KITTI_0, keypoints_KITTI);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	int n_keypoints = keypoints_KITTI.size();
	float mean_runtime = elapsed.count() / n_keypoints;
	std::cout << std::endl
		<< "Detector type: " << typeid(type).name() << std::endl
		<< "Number of keypoints: " << n_keypoints << '\n'
		<< "Elapsed time: " << elapsed.count() << " s\n"
		<< "Mean computation time per keypoint: " << mean_runtime << " s\n";
	return mean_runtime;
}

// Function 'mean_time_descriptors'
template<class type> float mean_time_descriptors(Ptr<type> &descriptor) {
	Ptr<FastFeatureDetector> detector_FAST = FastFeatureDetector::create();
	std::vector<KeyPoint> keypoints_KITTI;
	Mat descriptors_KITTI;
	detector_FAST->detect(img_KITTI_0, keypoints_KITTI);
	auto start = std::chrono::high_resolution_clock::now();
	descriptor->compute(img_KITTI_0, keypoints_KITTI, descriptors_KITTI);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	int n_keypoints = keypoints_KITTI.size();
	float mean_runtime = elapsed.count() / n_keypoints;
	std::cout << std::endl
		<< "Descriptor type: " << typeid(type).name() << std::endl
		<< "Number of keypoints: " << n_keypoints << '\n'
		<< "Elapsed time: " << elapsed.count() << " s\n"
		<< "Mean computation time per keypoint: " << mean_runtime << " s\n";
	return mean_runtime;
}

// Function 'mean_time_matchers'
template<class type> float mean_time_matchers(Ptr<type> &descriptor) {
	Ptr<FastFeatureDetector> detector_FAST = FastFeatureDetector::create();
	std::vector<KeyPoint> keypoints_KITTI_0, keypoints_KITTI_1;
	Mat  descriptors_KITTI_0, descriptors_KITTI_1;
	BFMatcher matcher;
	std::vector< DMatch > matches_KITTI;
	detector_FAST->detect(img_KITTI_0, keypoints_KITTI_0);
	detector_FAST->detect(img_KITTI_1, keypoints_KITTI_1);
	descriptor->compute(img_KITTI_0, keypoints_KITTI_0, descriptors_KITTI_0);
	descriptor->compute(img_KITTI_1, keypoints_KITTI_1, descriptors_KITTI_1);
	auto start = std::chrono::high_resolution_clock::now();
	matcher.match(descriptors_KITTI_0, descriptors_KITTI_1, matches_KITTI);
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	int n_keypoints = keypoints_KITTI_0.size() + keypoints_KITTI_1.size();
	int n_matches = matches_KITTI.size();
	float mean_runtime = elapsed.count() / n_matches;
	std::cout << std::endl
		<< "Descriptor type: " << typeid(type).name() << std::endl
		<< "Number of keypoints: " << n_keypoints << '\n'
		<< "Elapsed time: " << elapsed.count() << " s\n"
		<< "Mean computation time per match: " << mean_runtime << " s\n";
	return mean_runtime;
}

// Function 'main'
int main(int argc, char** argv)
{
	// Check for image loading errors
	if (!img_KITTI_0.data || !img_KITTI_1.data)
	{
		std::cout << " --(!) Error reading images " << std::endl; return -1;
	}

	// Initialize feature detectors:
	Ptr<FastFeatureDetector> detector_FAST = FastFeatureDetector::create();

	// Initialize feature descriptors:
	//Ptr<BRISK> descriptor_BRISK = BRISK::create();
	Ptr<ORB> descriptor_ORB = ORB::create();
	// Compute elapsed and mean runtime per keypoint detection
	int iterations = 10;
	std::ofstream results_detectors, results_descriptors, results_matchers;

	results_detectors.open("results_detectors.txt");
	for (int i = 0; i<iterations; i++) {
		results_detectors << mean_time_detectors<FastFeatureDetector>(detector_FAST) << '\t';
	
	}
	results_detectors.close();

	// Compute elapsed and mean runtime per keypoint description
	results_descriptors.open("results_descriptors.txt");
	for (int i = 0; i<iterations; i++) {
		results_descriptors /*<< mean_time_descriptors<SIFT>(descriptor_SIFT) << '\t'
							<< mean_time_descriptors<SURF>(descriptor_SURF) << '\t'*/
							//<< mean_time_descriptors<BriefDescriptorExtractor>(descriptor_BRIEF) << '\t'
			//<< mean_time_descriptors<BRISK>(descriptor_BRISK) << '\t'
			<< mean_time_descriptors<ORB>(descriptor_ORB) << '\t';
		/*<< mean_time_descriptors<FREAK>(descriptor_FREAK) << '\n';*/
	}
	results_descriptors.close();

	// Compute elapsed and mean runtime per keypoint match
	results_matchers.open("results_matchers.txt");
	for (int i = 0; i < iterations; i++) {
		results_matchers /*<< mean_time_matchers<SIFT>(descriptor_SIFT) << '\t'
						 << mean_time_matchers<SURF>(descriptor_SURF) << '\t'*/
						 //<< mean_time_matchers<BriefDescriptorExtractor>(descriptor_BRIEF) << '\t'
			//<< mean_time_matchers<BRISK>(descriptor_BRISK) << '\t'
			<< mean_time_matchers<ORB>(descriptor_ORB) << '\t';
		/*<< mean_time_matchers<FREAK>(descriptor_FREAK) << '\n';*/
	}
	results_matchers.close();

	return 0;
}