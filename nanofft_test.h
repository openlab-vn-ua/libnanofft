#pragma once
#ifndef NANOFFT_TEST_H
#define NANOFFT_TEST_H

extern bool testNanoFFT256(int directFFTCount = 1); // speed test
extern bool testNanoFFT4096(int directFFTCount = 1); // speed test
extern bool testNanoFFTAll(bool isLog = false);

#endif
