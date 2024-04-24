#pragma once
#ifndef NANOFFT_H
#define NANOFFT_H

// Fast Fourier transform proccessor
// --------------------------------------------------------------------------
// Proccess Direct and inverce fourier transform operations

// Code hardly based on free implemenation from Wikibooks:
// https://ru.wikibooks.org/wiki/%D0%A0%D0%B5%D0%B0%D0%BB%D0%B8%D0%B7%D0%B0%D1%86%D0%B8%D0%B8_%D0%B0%D0%BB%D0%B3%D0%BE%D1%80%D0%B8%D1%82%D0%BC%D0%BE%D0%B2/%D0%91%D1%8B%D1%81%D1%82%D1%80%D0%BE%D0%B5_%D0%BF%D1%80%D0%B5%D0%BE%D0%B1%D1%80%D0%B0%D0%B7%D0%BE%D0%B2%D0%B0%D0%BD%D0%B8%D0%B5_%D0%A4%D1%83%D1%80%D1%8C%D0%B5
// Ported to C++ without any serious changes of orginal code, may run even on Arduino platform
// Speed ~15% lower then KISSFFT library in simple float mode
// The main reason for this code is to have simple lightweight single .h file FFT implementaion

class NanoFFT
{
    // Fast Fourier Transform
    // PARAMETERS:  
    //    float *Rdat        [in, out] - Real part of Input and Output Data (Signal or Spectrum)
    //    float *Idat        [in, out] - Imaginary part of Input and Output Data (Signal or Spectrum)
    //    int    len         [in]      - Input and Output Data length (Number of samples in arrays)
    //    bool   isDirectFFT [in]      - true = Direct FFT (Signal to Spectrum), false = inverce FFT (Spectrum to Signal)
    //
    // RETURN VALUE: false on parameter error
    //
    // NOTE: In this algorithm N can be only power of 2 and in range:
    //       N    = 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384;
    //       LogN = 2, 3,  4,  5,  6,   7,   8,   9,   10,   11,   12,   13,    14;
    //

    // FFT assumes that samples array contains one period of periodic sample signal, see details:
    // https://gtest.com.ua/stati/osnovy-analiza-spektra-v-realnom-masshtabe-vremeni-chast-3.html

    // In short:
    // After Direct FFT:
    // DC part * COUNT will be at Re[0] (Im[0] will be 0)
    // Full 1.0 amp of [cos/sin] will be in [Re/Im] * COUNT * 0.5
    // Half of array will be progressively increasing freq cos/sin part
    // Re[1]/Im[1] will hold frequency that have one full cos/sin cycle over COUNT samples (1st harmonic)
    // Re[2]/Im[2] will hold frequency that have two full cos/sin cycle over COUNT samples (2nd harmonic)
    // After half of the array, the freq will be inverted (N-1 will be first harmonic), N-2 will 1/2 harmonic etc (?)

    public: 
    template<typename TINT> 
    static constexpr bool NUMBER_IS_2_POW_K(TINT x) { return ((!((x) & ((x)-1))) && ((x) > 1)); } // x is pow(2, k), k=1,2, ...

    public:
    template<typename TINT>
    static constexpr TINT GET_POWER_OF_2(int exp) { return (((TINT)1) << (exp)); }

    protected: static constexpr int UINT_LOG_2_INT(unsigned int val)
    {
        for (int i = 0; i < sizeof(val) * 8; i++)
        {
            if (val == GET_POWER_OF_2<decltype(val)>(i)) { return i; }
        }

        return -1;
    }

    public: static constexpr size_t FFT_BLOCK_MIN_SIZE = 4;
    public: static constexpr size_t FFT_BLOCK_MAX_SIZE = 16384;

    public: static bool FFT(float* Rdat, float* Idat, size_t len, bool isDirectFFT = true)
    {
        int N = len; // just to speedup, use shortest type possible, as N will fit into int
        // parameters error check:
        if ((Rdat == nullptr) || (Idat == nullptr))            return false;
        if ((N > 16384) || (N < 4))                            return false;
        if (!NUMBER_IS_2_POW_K(N))                             return false;
        int LogN = UINT_LOG_2_INT(N);
        if ((LogN < 2) || (LogN > 14))                         return false;

        int    i, j, n, k, io, ie, in, nn;
        float  ru, iu, rtp, itp, rtq, itq, rw, iw, sr;

        static const float Rcoef[14] =
        { 
           -1.0000000000000000F,  0.0000000000000000F,  0.7071067811865475F,
            0.9238795325112867F,  0.9807852804032304F,  0.9951847266721969F,
            0.9987954562051724F,  0.9996988186962042F,  0.9999247018391445F,
            0.9999811752826011F,  0.9999952938095761F,  0.9999988234517018F,
            0.9999997058628822F,  0.9999999264657178F
        };
        static const float Icoef[14] =
        { 
            0.0000000000000000F, -1.0000000000000000F, -0.7071067811865474F,
           -0.3826834323650897F, -0.1950903220161282F, -0.0980171403295606F,
           -0.0490676743274180F, -0.0245412285229122F, -0.0122715382857199F,
           -0.0061358846491544F, -0.0030679567629659F, -0.0015339801862847F,
           -0.0007669903187427F, -0.0003834951875714F
        };

        nn = N >> 1;
        ie = N;
        for (n = 1; n <= LogN; n++)
        {
            rw = Rcoef[LogN - n];
            iw = Icoef[LogN - n];
            if (!isDirectFFT) { iw = -iw; } // INVERSE
            in = ie >> 1;
            ru = 1.0F;
            iu = 0.0F;
            for (j = 0; j < in; j++)
            {
                for (i = j; i < N; i += ie)
                {
                    io = i + in;
                    rtp = Rdat[i] + Rdat[io];
                    itp = Idat[i] + Idat[io];
                    rtq = Rdat[i] - Rdat[io];
                    itq = Idat[i] - Idat[io];
                    Rdat[io] = rtq * ru - itq * iu;
                    Idat[io] = itq * ru + rtq * iu;
                    Rdat[i] = rtp;
                    Idat[i] = itp;
                }

                sr = ru;
                ru = ru * rw - iu * iw;
                iu = iu * rw + sr * iw;
            }

            ie >>= 1;
        }

        for (j = i = 1; i < N; i++)
        {
            if (i < j)
            {
                io = i - 1;
                in = j - 1;
                rtp = Rdat[in];
                itp = Idat[in];
                Rdat[in] = Rdat[io];
                Idat[in] = Idat[io];
                Rdat[io] = rtp;
                Idat[io] = itp;
            }

            k = nn;

            while (k < j)
            {
                j = j - k;
                k >>= 1;
            }

            j = j + k;
        }

        if (isDirectFFT) return true;

        rw = 1.0F / N;

        for (i = 0; i < N; i++)
        {
            Rdat[i] *= rw;
            Idat[i] *= rw;
        }

        return true;
    }
};

#endif
