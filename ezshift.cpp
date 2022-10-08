/*  Copyright (c) 2013, Robert Wang, email: robertwgh (at) gmail.com
   All rights reserved. https://github.com/robertwgh/ezSIFT
   Some algorithms used in this code referred to:
   1. OpenCV: http://opencv.org/
   2. VLFeat: http://www.vlfeat.org/
   The SIFT algorithm was developed by David Lowe. More information can be found
   from: David G. Lowe, "Distinctive image features from scale-invariant
   keypoints," International Journal of Computer Vision, 60, 2 (2004), pp.
   91-110.
   Pay attention that the SIFT algorithm is patented. It is your responsibility
   to use the code in a legal way. Patent information: Method and apparatus for
   identifying scale invariant features in an image and use of same for locating
   an object in an image  David G. Lowe, US Patent 6,711,293 (March 23, 2004).
   Provisional application filed March 8, 1999. Asignee: The University of
   British Columbia.
   Revision history:
      September 15th, 2013: initial version.
      July 8th, 2014: fixed a bug in sample_2x in image.h. The bug only happened
   for image with odd width or height. July 2nd, 2018: code refactor.
*/

#include "ezsift.h"
#include "common.h"
#include "image.h"
#include "timer.h"
#include "vvector.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <list>

namespace ezsift {

// Init sift parameters
void init_sift_parameters(bool doubleFirstOctave, float contrast_threshold,
                          float edge_threshold, float match_NDDR_threshold)
{
    SIFT_IMG_DBL = doubleFirstOctave;
    SIFT_CONTR_THR = contrast_threshold;
    SIFT_CURV_THR = edge_threshold;
    SIFT_MATCH_NNDR_THR = match_NDDR_threshold;
}

// Set up first Octave
// doubleFirstOctave = true, firstOcative=-1;
// doubleFirstOctave = false, firstOcative=0;
void double_original_image(bool doubleFirstOctave)
{
    SIFT_IMG_DBL = doubleFirstOctave;
    return;
}

// Compute octaves to build Gaussian Pyramid.
int build_octaves(const Image<unsigned char> &image,
                  std::vector<Image<unsigned char>> &octaves, int firstOctave,
                  int nOctaves)
{
    int w = image.w;
    int h = image.h;
    if (firstOctave == -1) {
        w = image.w * 2;
        h = image.h * 2;
    }

    for (int i = 0; i < nOctaves; i++) {
        if (i == 0 && firstOctave == -1) {
            octaves[i] = image.upsample_2x();
        }
        else if ((i == 0 && firstOctave != -1) ||
                 (i == 1 && firstOctave == -1)) {
            octaves[i] = image;
        }
        else {
            octaves[i] = octaves[(i - 1)].downsample_2x();
        }
        w = w / 2;
        h = h / 2;
    }
    return 0;
}

// Improved Gaussian Blurring Function
int gaussian_blur(const Image<float> &in_image, Image<float> &out_image,
                  std::vector<float> coef1d)
{
    int w = in_image.w;
    int h = in_image.h;
    int gR = static_cast<int>(coef1d.size()) / 2;

    Image<float> img_t(h, w);
    row_filter_transpose(in_image.data, img_t.data, w, h, &coef1d[0], gR);
    row_filter_transpose(img_t.data, out_image.data, h, w, &coef1d[0], gR);

    return 0;
}

// Apply Gaussian row filter to image, then transpose the image.
int row_filter_transpose(float *src, float *dst, int w, int h, float *coef1d,
                         int gR)
{
    float *row_buf = new float[w + gR * 2];
    float *row_start;
    int elemSize = sizeof(float);

    float *srcData = src;
    float *dstData = dst + w * h - 1;
    float partialSum = 0.0f;
    float *coef = coef1d;
    float *prow;

    float firstData, lastData;
    for (int r = 0; r < h; r++) {
        row_start = srcData + r * w;
        memcpy(row_buf + gR, row_start, elemSize * w);
        firstData = *(row_start);
        lastData = *(row_start + w - 1);
        for (int i = 0; i < gR; i++) {
            row_buf[i] = firstData;
            row_buf[i + w + gR] = lastData;
        }

        prow = row_buf;
        dstData = dstData - w * h + 1;
        for (int c = 0; c < w; c++) {
            partialSum = 0.0f;
            coef = coef1d;

            for (int i = -gR; i <= gR; i++) {
                partialSum += (*coef++) * (*prow++);
            }

            prow -= 2 * gR;
            *dstData = partialSum;
            dstData += h;
        }
    }
    delete[] row_buf;
    row_buf = nullptr;

    return 0;
}

// Build Gaussian pyramid using recursive method.
// The first layer is downsampled from last octave, layer=nLayers.
// All n-th layer is Gaussian blur from (n-1)-th layer.
int build_gaussian_pyramid(std::vector<Image<unsigned char>> &octaves,
                           std::vector<Image<float>> &gpyr, int nOctaves,
                           int nGpyrLayers)
{
    int nLayers = nGpyrLayers - 3;
    std::vector<std::vector<float>> gaussian_coefs =
        compute_gaussian_coefs(nOctaves, nGpyrLayers);

    int w, h;
    for (int i = 0; i < nOctaves; i++) {
        w = octaves[i].w;
        h = octaves[i].h;
        for (int j = 0; j < nGpyrLayers; j++) {
            if (i == 0 && j == 0) {
                gpyr[0].init(w, h);
                gaussian_blur(octaves[0].to_float(), gpyr[0],
                              gaussian_coefs[j]);
            }
            else if (i > 0 && j == 0) {
                gpyr[i * nGpyrLayers] =
                    gpyr[(i - 1) * nGpyrLayers + nLayers].downsample_2x();
            }
            else {
                gpyr[i * nGpyrLayers + j].init(w, h);
                gaussian_blur(gpyr[i * nGpyrLayers + j - 1],
                              gpyr[i * nGpyrLayers + j], gaussian_coefs[j]);
            }
        }
    }
    // Release octaves memory.
    octaves.clear();
    return 0;
}

// For build_gaussian_pyramid()
std::vector<std::vector<float>> compute_gaussian_coefs(int nOctaves,
                                                       int nGpyrLayers)
{
    // Compute all sigmas for different layers
    int nLayers = nGpyrLayers - 3;
    float sigma, sigma_pre;
    float sigma0 = SIFT_SIGMA;
    float k = powf(2.0f, 1.0f / nLayers);

    std::vector<float> sig(nGpyrLayers);
    sigma_pre = SIFT_IMG_DBL ? 2.0f * SIFT_INIT_SIGMA : SIFT_INIT_SIGMA;
    sig[0] = sqrtf(sigma0 * sigma0 - sigma_pre * sigma_pre);
    for (int i = 1; i < nGpyrLayers; i++) {
        sigma_pre = powf(k, (float)(i - 1)) * sigma0;
        sigma = sigma_pre * k;
        sig[i] = sqrtf(sigma * sigma - sigma_pre * sigma_pre);
    }

    std::vector<std::vector<float>> gaussian_coefs(nGpyrLayers);
    for (int i = 0; i < nGpyrLayers; i++) {
        // Compute Gaussian filter coefficients
        float factor = SIFT_GAUSSIAN_FILTER_RADIUS;
        int gR = (sig[i] * factor > 1.0f) ? (int)ceilf(sig[i] * factor) : 1;
        int gW = gR * 2 + 1;

        gaussian_coefs[i].resize(gW);
        float accu = 0.0f;
        float tmp;
        for (int j = 0; j < gW; j++) {
            tmp = (float)((j - gR) / sig[i]);
            gaussian_coefs[i][j] = expf(tmp * tmp * -0.5f) * (1 + j / 1000.0f);
            accu += gaussian_coefs[i][j];
        }
        for (int j = 0; j < gW; j++) {
            gaussian_coefs[i][j] = gaussian_coefs[i][j] / accu;
        } // End compute Gaussian filter coefs
    }
    return gaussian_coefs;
}

// Build difference of Gaussian pyramids.
int build_dog_pyr(std::vector<Image<float>> &gpyr,
                  std::vector<Image<float>> &dogPyr, int nOctaves,
                  int nDogLayers)
{
    int nGpyrLayers = nDogLayers + 1;

    int w, h;
    float *srcData1; // always data2-data1.
   
    for (int i = 0; i < nOctaves; i++) {
        int row_start = i * nGpyrLayers;
        w = gpyr[row_start].w;
        h = gpyr[row_start].h;

        for (int j = 0; j < nDogLayers; j++) {
            dogPyr[i * nDogLayers + j].init(w, h);
            dstData = dogPyr[i * nDogLayers + j].data;

            srcData1 = gpyr[row_start + j].data;
            srcData2 = gpyr[row_start + j + 1].data;

            index = 0;
            while (index++ < w * h)
                *(dstData++) = *(srcData2++) - *(srcData1++);
        }
    }

    return 0;
}
