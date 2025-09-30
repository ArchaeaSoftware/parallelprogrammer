/*
 * Copyright (C) 2025 by Nicholas Wilt.
 *
 * Redistribution and use in source and binary forms, with or without modification, 
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this 
 *    list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice, this
 *    list of conditions and the following disclaimer in the documentation and/or 
 *    other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may
 *    be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <string.h>
#include <strings.h>

#include <iostream>
#include <utility>
#include <algorithm>
#include <vector>
#include <limits>
#include <array>
#include <map>
#include <chrono>

#include <x86intrin.h>

template<typename T, std::size_t k>
inline void
SiftDown( std::array<T,k>& x, size_t i=0 )
{
    size_t c;
    while ( (c = i+i+1) < k ) {
        c += (c+1)<k && (x[c+1]<x[c]);
        //if ( (c+1)<k && (x[c+1]<x[c]) )
        //    c += 1;
        if ( x[i] < x[c] ) break;
        std::swap( x[i], x[c] );
        i = c;
    }
}

template<typename T, std::size_t k>
inline void
SiftUp( std::array<T,k>& x, size_t i )
{
    while ( i ) {
        size_t p = (i-1)>>1;
        if ( x[p] < x[i] )
            break;
        std::swap( x[i], x[p] );
        i = p;
    }
}


template<uint32_t k>
int32_t
kthLargest_heap( const std::vector<int32_t>& v )
{
    std::array<int32_t,k> heap;

    for ( size_t i = 0; i < k; i++ ) {
        heap[i] = v[i];
        SiftUp( heap, i );
    }
    int32_t minMax = heap[0];
    for ( size_t i = k; i < v.size(); i++ ) {
        int32_t x = v[i];

        if ( minMax < x ) {
            heap[0] = x;
            SiftDown( heap );
            minMax = heap[0];
        }
    }
    return heap[0];
}

template<uint32_t k>
int32_t
kthLargest_sort( const std::vector<int32_t>& v )
{
    // Local array we keep sorted in increasing order
    // The first element is smallest, so any incoming element
    // that is larger must be replace that element, then get
    // moved into position
    std::array<int32_t,k> maxSoFar;

    for ( size_t i = 0; i < k; i++ ) {
        maxSoFar[i] = v[i];
    }
    std::sort( maxSoFar.begin(), maxSoFar.end() );
    int32_t minMax = maxSoFar[0];
    for ( size_t i = k; i < v.size(); i++ ) {
        int32_t x = v[i];

        if ( minMax < x ) {
            size_t j;
            for ( j = 1; j < k && maxSoFar[j] < x; j++ ) {
                maxSoFar[j-1] = maxSoFar[j];
            }
            maxSoFar[j-1] = x;
            minMax = maxSoFar[0];
        }
    }
    return maxSoFar[0];
}


#if 0
template<bool fancySwap>
int32_t
thirdLargest( const std::vector<int32_t>& in )
{
    std::array<int32_t,3> v;

    auto swap_if = []( int32_t& x, int32_t& y, bool predicate ) -> void {
        if ( fancySwap ) {
            int32_t mask = -int32_t(predicate);
            x ^= y;
            y ^= (x&mask);
            x ^= y;
        }
        else {
            if ( predicate ) std::swap( x, y );
        }
    };

    auto insertNewMax = [&v, swap_if]( int32_t newValue ) -> void {
        v[0] = newValue;
        swap_if( v[0], v[1], v[0]>v[1] );
        swap_if( v[1], v[2], v[1]>v[2] );
    };
    {
        v[0] = in[0];
        v[1] = in[1];
        v[2] = in[2];

        swap_if( v[0], v[2], v[0]>v[2] );
        swap_if( v[0], v[1], v[0]>v[1] );
        swap_if( v[1], v[2], v[1]>v[2] );
    }

    for ( size_t i = 3; i < in.size(); ++i ) {
        if ( in[i] > v[0] ) {
            insertNewMax( in[i] );
        }
    }
    return v[0];
}
#endif

template<bool fancySwap>
int32_t
thirdLargest( const std::vector<int32_t>& v )
{
    int32_t v0, v1, v2;

    auto swap_if = []( int32_t& x, int32_t& y, bool predicate ) -> void {
        if ( fancySwap ) {
            int32_t mask = -int32_t(predicate);
            x ^= y;
            y ^= (x&mask);
            x ^= y;
        }
        else {
            if ( predicate ) std::swap( x, y );
        }
    };

    auto insertNewMax = [&v0, &v1, &v2, swap_if]( int32_t v3 ) -> void {
        v0 = v3;
        swap_if( v0, v1, v0>v1 );
        swap_if( v1, v2, v1>v2 );
    };
    {
        v0 = v[0];
        v1 = v[1];
        v2 = v[2];

        swap_if( v0, v2, v0>v2 );
        swap_if( v0, v1, v0>v1 );
        swap_if( v1, v2, v1>v2 );
    }

    for ( size_t i = 3; i < v.size(); ++i ) {
        if ( v[i] > v0 ) {
            insertNewMax( v[i] );
        }
    }
    return v0;
}

inline void insertNewMax_x4 ( __m128i& v_maxSoFar, __m128i v_x ) {
    __m128i v_cmplt = _mm_cmplt_epi32( v_maxSoFar, v_x );
    __m128i v_insertion_mask = _mm_xor_si128( v_cmplt, _mm_srli_si128( v_cmplt, 4 ) );
    __m128i v_shift1 = _mm_srli_si128( v_maxSoFar, 4 );

    v_maxSoFar = _mm_blendv_epi8( v_maxSoFar, v_shift1, v_cmplt );
    v_maxSoFar = _mm_blendv_epi8( v_maxSoFar, v_x, v_insertion_mask );
}

inline void insertNewMax( __m128i& v_maxSoFar, int32_t x ) {
    insertNewMax_x4( v_maxSoFar, _mm_set1_epi32( x ) );
}

inline void insertNewMax_x8 ( __m256i& v_maxSoFar, __m256i v_x ) {
    auto _mm256_srli_si256_4 = []( __m256i v_x ) -> __m256i {
        return _mm256_alignr_epi8(
            _mm256_permute2x128_si256(v_x, v_x, _MM_SHUFFLE(2, 0, 0, 1)), v_x, 4);
    };

    __m256i v_cmplt = _mm256_cmpgt_epi32( v_x, v_maxSoFar );
    __m256i v_insertion_mask = _mm256_xor_si256( v_cmplt, _mm256_srli_si256_4( v_cmplt ) );
    __m256i v_shift1 = _mm256_srli_si256_4( v_maxSoFar );
    v_maxSoFar = _mm256_blendv_epi8( v_maxSoFar, v_shift1, v_cmplt );
    v_maxSoFar = _mm256_blendv_epi8( v_maxSoFar, v_x, v_insertion_mask );
}


int32_t
thirdLargest_avx2( const std::vector<int32_t>& v )
{
    int32_t minMax;
    __m128i v_maxSoFar;

    {
        v_maxSoFar = _mm_set_epi32( std::numeric_limits<int32_t>::max(), 
                                    std::numeric_limits<int32_t>::min(),
                                    std::numeric_limits<int32_t>::min(),
                                    std::numeric_limits<int32_t>::min() );
        insertNewMax( v_maxSoFar, v[0] );
        insertNewMax( v_maxSoFar, v[1] );
        insertNewMax( v_maxSoFar, v[2] );
    }

    minMax = _mm_cvtsi128_si32( v_maxSoFar );
    size_t i = 3;
    while ( i < v.size() ) {
        if ( v[i] > minMax ) {
            insertNewMax( v_maxSoFar, v[i] );
            minMax = _mm_cvtsi128_si32( v_maxSoFar );
        }
        ++i;
    }

    return minMax;
}

int32_t
thirdLargest_avx2_x8( const std::vector<int32_t>& v )
{
    int32_t minMax;
    __m128i v_maxSoFar;

    {
        v_maxSoFar = _mm_set_epi32( std::numeric_limits<int32_t>::max(), 
                                    std::numeric_limits<int32_t>::min(),
                                    std::numeric_limits<int32_t>::min(),
                                    std::numeric_limits<int32_t>::min() );
        insertNewMax( v_maxSoFar, v[0] );
        insertNewMax( v_maxSoFar, v[1] );
        insertNewMax( v_maxSoFar, v[2] );
    }

    minMax = _mm_cvtsi128_si32( v_maxSoFar );
    size_t i = 3;
    __m256i v_minMax = _mm256_broadcastd_epi32( v_maxSoFar );
    size_t Nleft = v.size() - 3;
    while ( Nleft >= 8 ) {
        __m256i v_next = _mm256_loadu_si256( (__m256i *) (v.data() + i) );
        __m256i v_cmplt = _mm256_cmpgt_epi32( v_next, v_minMax );
        int mask_lt = _mm256_movemask_ps( _mm256_castsi256_ps( v_cmplt ) );

        auto insert_lane = [&v_maxSoFar, v_next]( int i ) {
            __m128i v_nexti = _mm256_castsi256_si128( _mm256_permutevar8x32_epi32( v_next, _mm256_broadcastd_epi32( _mm_cvtsi32_si128( i ) ) ) );
            insertNewMax_x4( v_maxSoFar, v_nexti );
        };

        if ( mask_lt ) {
            insert_lane( 0 );
            insert_lane( 1 );
            insert_lane( 2 );
            insert_lane( 3 );
            insert_lane( 4 );
            insert_lane( 5 );
            insert_lane( 6 );
            insert_lane( 7 );
        }
        v_minMax = _mm256_broadcastd_epi32( v_maxSoFar );
        Nleft -= 8;
        i += 8;
    }
    minMax = _mm_cvtsi128_si32( v_maxSoFar );
    while ( i < v.size() ) {
        if ( v[i] > minMax ) {
            insertNewMax( v_maxSoFar, v[i] );
            minMax = _mm_cvtsi128_si32( v_maxSoFar );
        }
        ++i;
    }

    return minMax;
}


int32_t
thirdLargest_avx2_x8_peel( const std::vector<int32_t>& v )
{
    int32_t minMax;
    __m128i v_maxSoFar;

    {
        v_maxSoFar = _mm_set_epi32( std::numeric_limits<int32_t>::max(), 
                                    std::numeric_limits<int32_t>::min(),
                                    std::numeric_limits<int32_t>::min(),
                                    std::numeric_limits<int32_t>::min() );
        insertNewMax( v_maxSoFar, v[0] );
        insertNewMax( v_maxSoFar, v[1] );
        insertNewMax( v_maxSoFar, v[2] );
    }

    minMax = _mm_cvtsi128_si32( v_maxSoFar );
    size_t i = 3;
    __m256i v_minMax = _mm256_broadcastd_epi32( v_maxSoFar );
    size_t Nleft = v.size() - 3;
    while ( Nleft >= 8 ) {
        __m256i v_next = _mm256_loadu_si256( (__m256i *) (v.data() + i) );
        __m256i v_cmplt = _mm256_cmpgt_epi32( v_next, v_minMax );
        int mask_lt = _mm256_movemask_ps( _mm256_castsi256_ps( v_cmplt ) );

        int mask_done = -1;
        while ( 0 != (mask_lt&mask_done) ) {
            int lsb = ffs( mask_done & mask_lt );
            __m128i v_newCandidate = _mm256_castsi256_si128( _mm256_permutevar8x32_epi32( v_next, _mm256_broadcastd_epi32( _mm_cvtsi32_si128( lsb-1) ) ) );
            insertNewMax_x4( v_maxSoFar, v_newCandidate );
            mask_done = ~((1<<lsb)-1);
        }
        v_minMax = _mm256_broadcastd_epi32( v_maxSoFar );
        Nleft -= 8;
        i += 8;
    }
    minMax = _mm_cvtsi128_si32( v_maxSoFar );
    while ( i < v.size() ) {
        if ( v[i] > minMax ) {
            insertNewMax( v_maxSoFar, v[i] );
            minMax = _mm_cvtsi128_si32( v_maxSoFar );
        }
        ++i;
    }

    return minMax;
}

template<uint32_t k>
int32_t
thirdLargest_avx2_x8_least( const std::vector<int32_t>& v )
{
    static_assert( k <= 8 );

    __m256i v_maxSoFar;
    {
        constexpr int32_t minInt = std::numeric_limits<int32_t>::min();
        constexpr int32_t maxInt = std::numeric_limits<int32_t>::max();
        int32_t v_init[8];
        for ( size_t i = 0; i < 8; i++ ) {
            v_init[i] = i<k ? minInt : maxInt;
        }
        v_maxSoFar = _mm256_loadu_si256( (__m256i *) &v_init );
        for ( size_t i = 0; i < k; i++ ) {
            insertNewMax_x8( v_maxSoFar, _mm256_set1_epi32( v[i] ) );
        }
    }

    int32_t minMax = _mm_cvtsi128_si32( _mm256_castsi256_si128( v_maxSoFar ) );
    __m256i v_minMax = _mm256_set1_epi32( minMax );
    size_t i = k;
    size_t Nleft = v.size() - k;
    while ( Nleft >= 8 ) {
        __m256i v_next = _mm256_loadu_si256( (__m256i *) (v.data() + i) );
        __m256i v_cmplt = _mm256_cmpgt_epi32( v_next, v_minMax );
        int mask_lt = _mm256_movemask_ps( _mm256_castsi256_ps( v_cmplt ) );

        int lsb = ffs( mask_lt );
        if ( lsb ) {
            __m256i v_newCandidate = _mm256_permutevar8x32_epi32( v_next, _mm256_broadcastd_epi32( _mm_cvtsi32_si128( lsb-1) ) );
            insertNewMax_x8( v_maxSoFar, v_newCandidate );
            minMax = _mm_cvtsi128_si32( _mm256_castsi256_si128( v_maxSoFar ) );
            v_minMax = _mm256_set1_epi32( minMax );

            Nleft -= lsb;
            i += lsb;
        }
        else {
            Nleft -= 8;
            i += 8;
        }
    }
    minMax = _mm_cvtsi128_si32( _mm256_castsi256_si128( v_maxSoFar ) );
    while ( i < v.size() ) {
        if ( v[i] > minMax ) {
            insertNewMax_x8( v_maxSoFar, _mm256_set1_epi32( v[i] ) );
            minMax = _mm_cvtsi128_si32( _mm256_castsi256_si128( v_maxSoFar ) );
        }
        ++i;
    }

    return minMax;
}

#if 0
template<uint32_t k>
int32_t
kthLargest_avx2_x8_least( const std::vector<int32_t>& v )
{
    static_assert( k < 8 );

    __m256i v_maxSoFar;
    {
        constexpr int32_t minInt = std::numeric_limits<int32_t>::min();
        constexpr int32_t maxInt = std::numeric_limits<int32_t>::max();
        int32_t v_init[8];
        for ( size_t i = 0; i < 8; i++ ) {
            v_init[i] = i<k ? maxInt : minInt;
        }
        v_maxSoFar = _mm256_loadu_si256( (__m256i *) &v_init );
        for ( size_t i = 0; i < k; i++ ) {
            insertNewMax_x8( v_maxSoFar, _mm256_set1_epi32( v[i] ) );
        }
    }

    int32_t minMax = _mm_cvtsi128_si32( _mm256_castsi256_si128( v_maxSoFar ) );
    __m256i v_minMax = _mm256_set1_epi32( minMax );
    size_t i = k;
    size_t Nleft = v.size() - k;
    while ( Nleft >= 8 ) {
        __m256i v_next = _mm256_loadu_si256( (__m256i *) (v.data() + i) );
        __m256i v_cmplt = _mm256_cmpgt_epi32( v_next, v_minMax );
        int mask_lt = _mm256_movemask_ps( _mm256_castsi256_ps( v_cmplt ) );

        int lsb = ffs( mask_lt );
        if ( lsb ) {
            __m256i v_newCandidate = _mm256_permutevar8x32_epi32( v_next, _mm256_broadcastd_epi32( _mm_cvtsi32_si128( lsb-1) ) );
            insertNewMax_x8( v_maxSoFar, v_newCandidate );
            minMax = _mm_cvtsi128_si32( _mm256_castsi256_si128( v_maxSoFar ) );
            v_minMax = _mm256_set1_epi32( minMax );

            Nleft -= lsb;
            i += lsb;
        }
        else {
            Nleft -= 8;
            i += 8;
        }
    }
    while ( i < v.size() ) {
        if ( v[i] > minMax ) {
            insertNewMax_x8( v_maxSoFar, _mm256_set1_epi32( v[i] ) );
            minMax = _mm_cvtsi128_si32( _mm256_castsi256_si128( v_maxSoFar ) );
        }
        ++i;
    }

    return minMax;
}
#endif

int
main( int argc, char *argv[] )
{
    bool report_k = false;
    bool report_csv = false;
    bool sortInputs = false;
    bool revSortInputs = false;
    bool report_third = false;

    if ( 1 == argc ) {
        std::cout << "Please specify at least one command line argument:" << std::endl;
        std::cout << "    --sort: sort input array (best runtime)" << std::endl;
        std::cout << "    --revsort: reverse-sort input array (worst runtime)" << std::endl;
        std::cout << "    --csv: comma-delimited output" << std::endl;
        std::cout << "    --report_third: report performance for 3rd-largest variants" << std::endl;
        std::cout << "    --report_k: report performance for kth-largest variants" << std::endl;
        exit( 0 );
    }

    for ( int i = 1; i < argc; i++ ) {
        if ( 0 == strcmp( argv[i], "--sort" ) ) {
            sortInputs = true;
        }
        if ( 0 == strcmp( argv[i], "--revsort" ) ) {
            revSortInputs = true;
        }
        if ( 0 == strcmp( argv[i], "--csv" ) ) {
            report_csv = true;
        }
        if ( 0 == strcmp( argv[i], "--report_k" ) ) {
            report_k = true;
        }
        if ( 0 == strcmp( argv[i], "--report_third" ) ) {
            report_third = true;
        }
    }
    if ( sortInputs && revSortInputs ) {
        std::cout << "--sort and --revsort are mutually exclusive" << std::endl;
        exit( 1 );
    }
    if ( ! (report_k || report_third) ) {
        std::cout << "Please specify either --report_k or --report_third" << std::endl;
        exit( 1 );
    }

    constexpr size_t nIters = 10;
    constexpr size_t N = 100ULL*1048576;
    std::vector<int32_t> v( N );
    std::vector<int32_t> vOriginal( N );

    for ( auto& i : v ) {
        i = rand();
    }
    vOriginal = v;
    std::vector<int32_t> sorted( v );
    std::sort( sorted.begin(), sorted.end() );

    // maps for timings of various sizes
    std::map<size_t, std::array<double, nIters>> et_nth;
    std::map<size_t, std::array<double, nIters>> et_sort;
    std::map<size_t, std::array<double, nIters>> et_heap;

    auto zeroArray = []( ) -> std::array<double, nIters> {
        std::array<double, nIters> ret;
        for ( auto& d : ret ) d = 0.0f;
        return ret;
    };

    for ( size_t i = 100; i <= 1000; i+=100 ) {
        et_nth[i] = zeroArray( );
        et_sort[i] = zeroArray( );
        et_heap[i] = zeroArray( );
    }
    for ( size_t i = 2000; i <= 10000; i+=1000 ) {
        et_nth[i] = zeroArray( );
        et_sort[i] = zeroArray( );
        et_heap[i] = zeroArray( );
    }

    for ( size_t i = 0; i < nIters; i++ ) {
        for ( auto& i : v ) {
            i = rand();
        }
        if ( sortInputs ) {
            std::sort( v.begin(), v.end() );
        }
        else if ( revSortInputs ) {
            std::sort( v.begin(), v.end(), std::greater<int32_t>() );
        }

        vOriginal = v;
        sorted = v;
        std::sort( sorted.begin(), sorted.end() );
        int32_t expectedValue = sorted[sorted.size()-3];

        v = vOriginal;
        std::nth_element( v.begin(), v.begin()+v.size()-3, v.end() );
        assert( expectedValue == v[v.size()-3] );

        auto timeTest = [&v, vOriginal]( auto pfn, int32_t expectedValue ) -> double {
            v = vOriginal;
            auto start = std::chrono::high_resolution_clock::now();
            int32_t thirdLargestElement = pfn( v );
            assert( thirdLargestElement == expectedValue );
            auto stop = std::chrono::high_resolution_clock::now();
            return std::chrono::duration_cast<std::chrono::nanoseconds>(stop-start).count() / 1e9;
        };

        double et_nth_element;
        std::cout << "nth_element" << ((report_csv)?", ": ": ");
        std::cout << (et_nth_element = timeTest( []( auto& v ) -> uint64_t { 
                                            std::nth_element( v.begin(), v.begin()+v.size()-3, v.end() );
                                            return v[v.size()-3]; }, expectedValue ));
        std::cout << ((report_csv) ? ", 1.0":" s") << std::endl;

        auto reportTime = [timeTest, et_nth_element, report_csv]( auto pfn, const char *s, int32_t expectedValue, double *pet = NULL ) -> void {
            std::cout << s;
            if ( report_csv ) {
                std::cout << ", ";
            } else {
                std::cout << ": ";
            }
            auto et = timeTest( pfn, expectedValue );
            if ( pet ) *pet = (double) et;
            if ( report_csv ) {
                std::cout << et << ", " << et_nth_element/et << std::endl;
            }
            else {
                std::cout << et << " s (" << et_nth_element / et << "x faster)" << std::endl;
            }
        };

        if ( report_third ) {
            reportTime( kthLargest_heap<3>, "kthLargest_heap", expectedValue );
            reportTime( kthLargest_sort<3>, "kthLargest_sort", expectedValue );
            reportTime( thirdLargest<false>, "3rdLargest_plainswap", expectedValue );
            reportTime( thirdLargest<false>, "3rdLargest_fancyswap", expectedValue );
            reportTime( thirdLargest_avx2, "3rdLargest_avx2", expectedValue );
            reportTime( thirdLargest_avx2_x8, "3rdLargest_avx2_x8", expectedValue );
            reportTime( thirdLargest_avx2_x8_peel, "3rdLargest_avx2_x8_peel", expectedValue );
            reportTime( thirdLargest_avx2_x8_least<3>, "3rdLargest_avx2_x8_least", expectedValue );
        }

#define REPORTTIME(k) \
        std::cout << "nth_element" << ((report_csv) ? ", ":": ") << \
                     (et_nth_element = timeTest( []( auto& v ) -> uint64_t { \
                                            std::nth_element( v.begin(), v.begin()+v.size()-k, v.end() ); \
                                            return v[v.size()-k]; }, sorted[sorted.size()-k] )) << " s" << std::endl; \
        et_nth[k][i] = et_nth_element; \
        reportTime( kthLargest_heap<k>, "kthLargest_heap (k==" #k ")", sorted[sorted.size()-k], &et_heap[k][i] ); \
        reportTime( kthLargest_sort<k>, "kthLargest_sort (k==" #k ")", sorted[sorted.size()-k], &et_sort[k][i] );

        if ( report_k ) {
            REPORTTIME(100)
            REPORTTIME(200)
            REPORTTIME(300)
            REPORTTIME(400)
            REPORTTIME(500)
            REPORTTIME(600)
            REPORTTIME(700)
            REPORTTIME(800)
            REPORTTIME(900)
            REPORTTIME(1000)

            REPORTTIME(2000)
            REPORTTIME(3000)
            REPORTTIME(4000)
            REPORTTIME(5000)
            REPORTTIME(6000)
            REPORTTIME(7000)
            REPORTTIME(8000)
            REPORTTIME(9000)
            REPORTTIME(10000)

            REPORTTIME(20000)
            REPORTTIME(30000)
            REPORTTIME(40000)
            REPORTTIME(50000)
            REPORTTIME(60000)
            REPORTTIME(70000)
            REPORTTIME(80000)
            REPORTTIME(90000)
            REPORTTIME(100000)
        }

    }

    if ( report_k ) {
        std::map< size_t, double > average_nth;
        std::map< size_t, double > average_sort;
        std::map< size_t, double > average_heap;

        // compute average elements per second
        auto average = [N]( const std::array<double, nIters>& v ) -> double {
            double sum = 0.0;
            for ( size_t i = 0; i < nIters; ++i ) {
                sum += v[i];
            }
            sum /= (double) nIters;
            return (double) N / sum;
        };

        for ( auto v : et_nth ) {
            average_nth[v.first] = average( v.second );
            average_sort[v.first] = average( et_sort[v.first] );
            average_heap[v.first] = average( et_heap[v.first] );
        }
        std::cout << "k, nth, sort, heap" << std::endl;
        for ( auto v : average_nth ) {
            std::cout << v.first << ", " << v.second << ", " << average_sort[v.first] << ", " << average_heap[v.first] << std::endl;
        }
    }
    return 0;
}

