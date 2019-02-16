# special course

## week-1

### Content

1. wikipedia on rotation representation: quaternion, euler angle, rotation matrix, axis-angle; 附件旋转表达
2. useful tools:
   * LaTex for paper writing;
   * Github for code version control and review;

### Assignment

1. writing a library for conversion among different kinds of representations of rotation, e.g. q2r transforms a quaternion (input) to the corresponding rotation matrix. Ideal package should include (but not limited to):
   1. q2r: quaternion to rotation matrix;
   2. r2q: rotation matrix to quaternion;
   3. e2r: Euler angles to rotation matrix using ZYX order;
   4. r2e: rotation matrix to Euler angles using ZYX order;
   5. q2e: quaternion to Euler angles;
   6. ...... 

2. create a Github account and submit your implementation to your first repository.

### Reading materials:

1. representation of 3D rigid motion;
2. wikipedia on Euler angles, quaternion, rotation matrix;

## week-2 & 3

### Content

1. Fundamentals of Computer Vision: camera model and calibration principles (perspective camera);
2. Correspondence problem;
3. image features: corner and blob features.

### Assignment 

1. implement the Harris corner detector;
2. implement the Fast feature detector;

### Reading materials:

1. Introduction to Autonomous Mobile Robots: 
   * 4.2.1-4.2.3, 
   * 4.3.1, 
   * 4.4 Feature Extraction
   * 4.5 Image Feature Extraction: Interest Point Detectors: **Harris**, **FAST**.
2. https://www.edwardrosten.com/work/fast.html
3. camera and image (相机与图像).
4. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_3_2_0_introduction_to_keypoint_features_annotated.pdf
5. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_3_2_1_corner_features.pdf
6. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_3_2_2_blob_features.pdf

## week-4

### Content

1. feature descriptor and feature matching;
2. strategy: brute-force, k-nearest-neighbor;

### Assignment

1. Introduction to Autonomous Mobile Robots: 
   - 4.2.5: correspondence problem
   - 4.3.3 

1. extract FAST features and describe them using BRIEF descriptor.
2. feature matching using brute-force strategy.

### Reading material

1. https://www.cs.ubc.ca/~lowe/525/papers/calonder_eccv10.pdf: brief.
2. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_4_0_from_keypoints_to_correspondences.pdf
3. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_4_1_feature_descriptors.pdf
4. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_4_2_feature_matching.pdf

## week-5 

### Content

1. robust estimation: RANSAC;

### Assignment

1. Use RANSAC for line and plane estimation;

### Reading material

1. Introduction to Autonomous Mobile Robots: 
   - 4.7.2.3
2. http://www.cse.psu.edu/~rtc12/CSE486/lecture15.pdf
3. https://en.wikipedia.org/wiki/Random_sample_consensus
4. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_3_3-robust-estimation-with-ransac.pdf

## week-6

### Content

1. Two view geometry: essential, fundamental, **Homography**;
2. DLT for homography estimation;

### Assignment

1. extract features, matching, use DLT for homography estimation;
2. use RANSAC to find inliers and outliers;
3. image warping;

### Reading material

1. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_4_3-estimating-homographies-from-feature-correspondences.pdf
2. https://www.uio.no/studier/emner/matnat/its/UNIK4690/v16/forelesninger/lecture_6_1_basic_epipolar-geometry.pdf

## week-7 & 8

### Content

1. misc: filter, IMU, attitude estimation, etc (**only week 7**).

### Assignment

1. panarama: image-stitching
   * try with SIFT and FAST+BRIEF as feature for matching without RANSAC for outlier removal;
   * try with RANSAC and use homography for outlier removal;
2. Spherical image stiching: same with previous but project points into a sphere.

### Reading material

**TBD**

## week-9&10&11 

### Assignment

1. benchmark of feature detectors, descriptors and matching strategy;
2. benchmark different kind of matching outlier removal strategy;

### Reading material

**TBD**

## week-11-13

report writing.