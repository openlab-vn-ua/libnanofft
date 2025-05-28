#include "nanofft_test.h"
#include "nanofft.h"

#include <string.h>
#include <math.h>

#if !defined(logRawPrintln)
#include <iostream>
#define logRawPrintln(x) std::cout << (x) << "\n";
#define logRawPrintF(fmt,...) printf(fmt, __VA_ARGS__)
#else
#endif

static bool areFloatEqual(float a, float b, float EPS = 0.001)
{
    if (a == b) { return true; }
    if ((a >= b) && ((a - b) < EPS)) { return true; }
    if ((a <= b) && ((b - a) < EPS)) { return true; }
    return false;
}

template<int COUNT = 16>
static bool testNanoFFTStep(bool isLog, int directFFTCount = 1, int inverseFFTCount = 1)
{
    #define logPFX "testNanoFFTStep:"

    bool isOk = true;

    bool isLogVerbose = false;

    if (directFFTCount <= 0) { return false; } // invalid params
    //static constexpr int COUNT = 16; // number of samples
    constexpr float EPS = COUNT <= 256 ? 0.001 : 0.05;

    // Source
    float srcR[COUNT];
    float srcI[COUNT];
    float amp = 1.0;
    float bp = 2 * 3.141592653589 / COUNT; // There will be COUNT samples per period of base freq period
    float p = bp; // actual period (mul to increase freq)
    float amp0 = amp;
    float amp1 = amp;
    float amp2 = amp;
    float amp3 = amp;

    int i;
    for (i = 0; i < COUNT; i++)
    {
        srcR[i] = 0.0; // Real part of signal      (cos)
        srcI[i] = 0.0; // Imaginary part of signal (sin)

        srcR[i] += amp0 * 0.5; // DC part is doubled
        srcR[i] += amp1 * cos(p * i * 1);  // Real part of signal
        srcR[i] += amp2 * sin(p * i * 2);  // Real part of signal double freq (even harmonis will be sin in this test)
        srcR[i] += amp3 * cos(p * i * 3);  // Real part of signal tripple freq
        srcI[i] = 0.0;
    }

    if (isLog && isLogVerbose)
    {
        if (isLog) { logRawPrintF(logPFX "Source[%d]:\n", (int)COUNT); }
        for (i = 0; i < COUNT; i++)
        {
            if (isLog) { logRawPrintF(logPFX "%10.6f  %10.6f  %10.6f\n", srcR[i], srcI[i], srcR[i] * srcR[i] + srcI[i] * srcI[i]); }
        }
    }

    // Spectral
    float fftR[COUNT];
    float fftI[COUNT];

    for (i = 0; i < directFFTCount; i++)
    {
        memcpy(fftR, srcR, sizeof(fftR));
        memcpy(fftI, srcI, sizeof(fftI));
        NanoFFT::FFT(fftR, fftI, COUNT, true); // Do direct FFT
    }

    if (isLog && isLogVerbose)
    {
        // Output real and imaginary parts of the FFT and spectral power
        if (isLog) { logRawPrintF(logPFX "FFT[%d]:\n", (int)COUNT); }
        for (i = 0; i < COUNT; i++)
        {
            if (isLog) { logRawPrintF(logPFX "%10.6f  %10.6f  %10.6f\n", fftR[i], fftI[i], fftR[i] * fftR[i] + fftI[i] * fftI[i]); }
        }
    }

    if (!areFloatEqual(fftR[0],  0.5 * amp0 * COUNT, EPS)) { isOk = false; }
    if (!areFloatEqual(fftI[0],  0.0, EPS))                { isOk = false; }

    if (!areFloatEqual(fftR[1],  0.5 * amp1 * COUNT, EPS)) { isOk = false; }
    if (!areFloatEqual(fftI[1],  0.0, EPS))                { isOk = false; }

    if (!areFloatEqual(fftR[2],  0.0, EPS))                { isOk = false; }
    if (!areFloatEqual(fftI[2], -0.5 * amp2 * COUNT, EPS)) { isOk = false; } // negative here

    if (!areFloatEqual(fftR[3],  0.5 * amp3 * COUNT, EPS)) { isOk = false; }
    if (!areFloatEqual(fftI[3],  0.0, EPS))                { isOk = false; }

    for (i = 4; i < COUNT - 3; i++)
    {
        if (!areFloatEqual(fftR[i], 0.0, EPS)) { isOk = false; }
        if (!areFloatEqual(fftI[i], 0.0, EPS)) { isOk = false; }
    }

    if (!areFloatEqual(fftR[COUNT - 1],  0.5 * amp1 * COUNT, EPS)) { isOk = false; }
    if (!areFloatEqual(fftI[COUNT - 1],  0.0, EPS))                { isOk = false; }

    if (!areFloatEqual(fftR[COUNT - 2],  0.0, EPS))                { isOk = false; }
    if (!areFloatEqual(fftI[COUNT - 2],  0.5 * amp2 * COUNT, EPS)) { isOk = false; } // positive here

    if (!areFloatEqual(fftR[COUNT - 3],  0.5 * amp3 * COUNT, EPS)) { isOk = false; }
    if (!areFloatEqual(fftI[COUNT - 3],  0.0, EPS))                { isOk = false; }

    // Restored
    float outR[COUNT];
    float outI[COUNT];

    if (inverseFFTCount > 0)
    {
        for (i = 0; i < inverseFFTCount; i++)
        {
            memcpy(outR, fftR, sizeof(outR));
            memcpy(outI, fftI, sizeof(outI));
            NanoFFT::FFT(outR, outI, COUNT, false); // Do inverse FFT
        }

        if (isLog && isLogVerbose)
        {
            if (isLog) { logRawPrintF(logPFX "Inverse[%d]:\n", (int)COUNT); }
            for (i = 0; i < COUNT; i++)
            {
                if (isLog) { logRawPrintF(logPFX "%10.6f  %10.6f  %10.6f |-| %10.6f  %10.6f\n", outR[i], outI[i], outR[i] * outR[i] + outI[i] * outI[i], outR[i]-srcR[i], outI[i] - srcI[i]); }
            }
        }

        for (i = 0; i < COUNT; i++)
        {
            if (!areFloatEqual(srcR[i], outR[i], EPS)) { isOk = false; }
            if (!areFloatEqual(srcI[i], outI[i], EPS)) { isOk = false; }
        }
    }

    if (isLog)
    {
        logRawPrintF(logPFX "[%d]:%s" "\n", (int)COUNT, (isOk ? "OK" : "FAIL"));
    }

    #undef logPFX

    return isOk;
}

bool testNanoFFT256(int directFFTCount)
{
    return testNanoFFTStep<256>(false, directFFTCount, 0);
}

bool testNanoFFT4096(int directFFTCount)
{
    return testNanoFFTStep<4096>(false, directFFTCount, 0);
}

bool testNanoFFTAll(bool isLog)
{
    bool isOk = true;
    if (!testNanoFFTStep<8>(isLog)) { isOk = false; }
    if (!testNanoFFTStep<16>(isLog)) { isOk = false; }
    if (!testNanoFFTStep<256>(isLog)) { isOk = false; }
  //if (!testNanoFFTStep<4096>(isLog)) { isOk = false; }
    return isOk;
}
