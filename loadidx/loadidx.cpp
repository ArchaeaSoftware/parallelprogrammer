#include <x86intrin.h>
#include <stdint.h>
#include <stdio.h>

#include <asmjit/asmjit.h>

using namespace asmjit;

typedef __m256i (*pfnGetIndices)(void);

// Signature of the generated function.
typedef __m256i (*Func)(void);

int
emitIdentityIndices( JitRuntime& rt, Func *pfn )
{

  // Holds code and relocation information during code generation.
  CodeHolder code;

  // Code holder must be initialized before it can be used. The simples way to initialize
  // it is to use 'Environment' from JIT runtime, which matches the target architecture,
  // operating system, ABI, and other important properties.
  code.init(rt.environment(), CpuFeatures::X86::kAVX2 );//rt.cpuFeatures());

  // Emitters can emit code to CodeHolder - let's create 'x86::Assembler', which can emit
  // either 32-bit (x86) or 64-bit (x86_64) code. The following line also attaches the
  // assembler to CodeHolder, which calls 'code.attach(&a)' implicitly.
  x86::Assembler a(&code);

  // Use the x86::Assembler to emit some code to .text section in CodeHolder:

  x86::Xmm m8, m32hi;
  x86::Gp indices8, four;
  x86::Xmm v_fours, v_zeros;

  indices8 = x86::eax;
  four = x86::edx;
  m8 = x86::xmm1;
  m32hi = x86::xmm2;
  v_fours = x86::xmm3;
  v_zeros = x86::xmm4;

  a.mov( indices8, 0x03020100 );
  a.xor_( four, four ); // xor is a C++ keyword!
  a.mov( x86::dl, 4 );
  a.xorpd( v_zeros, v_zeros );
  a.movd( m8, indices8 );
#if NO_AVX512
  a.movd( v_fours, four ); // see if we can get rid of this and go straight from GPR
  a.vpbroadcastd( v_fours, v_fours );
#else
  a.vpbroadcastd( v_fours, four );
#endif
  a.vpunpcklbw( m8, m8, v_zeros );
  a.vpunpcklwd( m8, m8, v_zeros );
  a.vpaddd( v_fours, v_fours, m8 );
  a.vinserti128( x86::ymm0, m8, v_fours, 1 );
  a.ret();

  // 'x86::Assembler' is no longer needed from here and can be destroyed or explicitly
  // detached via 'code.detach(&a)' - which detaches an attached emitter from code holder.

  // Now add the generated code to JitRuntime via JitRuntime::add(). This function would
  // copy the code from CodeHolder into memory with executable permission and relocate it.
  Error err = rt.add(pfn, &code);

  // It's always a good idea to handle errors, especially those returned from the Runtime.
  if (err) {
    printf("AsmJit failed: %s\n", DebugUtils::errorAsString(err));
    return 1;
  }

  // CodeHolder is no longer needed from here and can be safely destroyed. The runtime now
  // holds the relocated function, which we have generated, and controls its lifetime. The
  // function will be freed with the runtime, so it's necessary to keep the runtime around.
  //
  // Use 'code.reset()' to explicitly free CodeHolder's content when necessary.

  // Execute the generated function and print the resulting '1', which it moves to 'eax'.
  __m256i result = (*pfn)();
  printf("%d\n", _mm_cvtsi128_si32( _mm256_castsi256_si128( result) ) );

  // All classes use RAII, all resources will be released before `main()` returns, the
  // generated function can be, however, released explicitly if you intend to reuse or
  // keep the runtime alive, which you should in a production-ready code.
  //rt.release(fn);

  return 0;
}

template<typename T>
void
print( const T& x )
{
    const char *s = (const char *) &x;
    printf( "%dB: 0x", (int) sizeof(T) );
    for ( int i = (int) (sizeof(T)-1); i >= 0; i-- ) {
        printf( "%02x", s[i] );
    }
    printf( "\n" );
}

__attribute((noinline))
__m256i get_indices1()
{
    return _mm256_set_epi32( 7, 6, 5, 4, 3, 2, 1, 0 );
}

// (removed get_indices9: dead, broken experiment superseded by get_indices10/11)


__attribute((noinline))
__m256i get_indices2()
{
    __m128i m8 = _mm_cvtsi32_si128( 0x03020100 );
    __m128i m16 = _mm_unpacklo_epi8( m8, _mm_setzero_si128() );
    __m128i m16hi = _mm_add_epi16( m16, _mm_set1_epi16( 4 ) );
    __m128i m32 = _mm_unpacklo_epi16( m16, _mm_setzero_si128() );
    __m128i m32hi = _mm_add_epi32( m32, _mm_set1_epi32( 4 ) );

    return _mm256_inserti128_si256( _mm256_castsi128_si256( m32 ), m32hi, 1 );
}


__attribute((noinline))
__m256i get_indices3()
{
    __m128i m8lo, m8hi;
    __m128i m32hi;
    uint32_t indices_low, indices_high;
    __m128i v_fours, v_zeros = _mm_setzero_si128();

    asm("movl $0x03020100, %0" : "=&r"(indices_low) );
    asm("movl $0x07060504, %0" : "=&r"(indices_high) );
    asm("movd %1, %0" : "=x"(m8lo) : "r"(indices_low) );
    asm("movd %1, %0" : "=x"(m8hi) : "r"(indices_high) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8lo) : "x"(v_zeros) );
    asm("vpunpcklwd %1, %0, %0" : "+x"(m8lo) : "x"(v_zeros) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8hi) : "x"(v_zeros) );
    asm("vpunpcklwd %1, %0, %0" : "+x"(m8hi) : "x"(v_zeros) );

    return _mm256_inserti128_si256( _mm256_castsi128_si256( m8lo ), m8hi, 1 );
}

__attribute((noinline))
__m256i get_indices4()
{
    __m128i m8;
    __m128i m32hi;
    uint32_t indices8, four;
    __m128i v_fours, v_zeros = _mm_setzero_si128();

    asm("movl $0x03020100, %0" : "=&r"(indices8) );
    asm("movl $4, %0" : "=&r"(four) );
    asm("movd %1, %0" : "=x"(m8) : "r"(indices8) );
    asm("movd %1, %0" : "=x"(v_fours) : "r"(four) );
    asm("vpbroadcastd %1, %0" : "=x"(v_fours) : "x"(v_fours) );
    //asm("vpbroadcastd %1, %0" : "=x"(v_fours) : "r"(four) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8) : "x"(v_zeros) );
    asm("vpunpcklwd %1, %0, %0" : "+x"(m8) : "x"(v_zeros) );
    asm("vpaddd %2, %1, %0" : "=x"(m32hi) : "x"(v_fours), "x"(m8) );

    return _mm256_inserti128_si256( _mm256_castsi128_si256( m8 ), m32hi, 1 );
}

__attribute((noinline))
__m256i get_indices5()
{
    __m128i m8lo, m8hi;
    __m128i m32hi;
    uint32_t indices_low, indices_high;
    __m128i v_fours, v_zeros = _mm_setzero_si128();

    asm("movl $0x03020100, %0" : "=&r"(indices_low) );
    asm("movl $0x07060504, %0" : "=&r"(indices_high) );
    asm("movd %1, %0" : "=x"(m8lo) : "r"(indices_low) );
    asm("movd %1, %0" : "=x"(m8hi) : "r"(indices_high) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8lo) : "x"(v_zeros) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8hi) : "x"(v_zeros) );

    __m256i m16 = _mm256_inserti128_si256( _mm256_castsi128_si256( m8lo ), m8hi, 1 );
    return _mm256_unpacklo_epi16( m16, _mm256_setzero_si256() );//, m16 );
}

__attribute((noinline))
__m256i get_indices6()
{
    __m128i m8lo, m8hi;
    __m128i m32hi;
    uint32_t indices_low, indices_high;
    __m128i v_fours, v_zeros = _mm_setzero_si128();
    static union {
        uint64_t indices = 0x0706050403020100ull;
        struct {
            uint32_t indiceslo;
            uint32_t indiceshi;
        };
    } c;

    asm("movd %1, %0" : "=x"(m8lo) : "m"(c.indiceslo) );
    asm("movd %1, %0" : "=x"(m8hi) : "m"(c.indiceshi) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8lo) : "x"(v_zeros) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8hi) : "x"(v_zeros) );

    __m256i m16 = _mm256_inserti128_si256( _mm256_castsi128_si256( m8lo ), m8hi, 1 );
    return _mm256_unpacklo_epi16( m16, _mm256_setzero_si256() );//, m16 );
}

__attribute((noinline))
__m256i get_indices7()
{
    __m128i m8;
    __m128i m32hi;
    uint32_t indices8;
    union { uint32_t four; uint8_t four0; };
    __m128i v_fours, v_zeros = _mm_setzero_si128();

    asm("xor %0,%0" : "=&r"(four), "=&r"(four) );
    asm("movl $0x03020100, %0" : "=&r"(indices8) );
    asm("movb $4, %0" : "=&r"(four0) );
    asm("movd %1, %0" : "=x"(m8) : "r"(indices8) );

    asm("movd %1, %0" : "=x"(v_fours) : "r"(four) );
    asm("vpbroadcastd %1, %0" : "=x"(v_fours) : "x"(v_fours) );

//    asm("vpbroadcastd %1, %0" : "=x"(v_fours) : "r"(four) );
    asm("vpunpcklbw %1, %0, %0" : "+x"(m8) : "x"(v_zeros) );
    asm("vpunpcklwd %1, %0, %0" : "+x"(m8) : "x"(v_zeros) );
    asm("vpaddd %2, %1, %0" : "=x"(m32hi) : "x"(v_fours), "x"(m8) );

    return _mm256_inserti128_si256( _mm256_castsi128_si256( m8 ), m32hi, 1 );
}


// The winner: build the 8-byte identity pattern in a GPR (immediate in the
// instruction stream), cross to XMM with vmovq, then let vpmovzxbd fan the 8
// bytes out to 8 dwords across the YMM in a single reg-to-reg op. No memory.
__attribute((noinline))
__m256i get_indices10()
{
    // _mm_cvtsi64_si128 of a constant -> movabs + vmovq (no memory load).
    __m128i m8 = _mm_cvtsi64_si128( 0x0706050403020100ll );
    // vpmovzxbd ymm, xmm : zero-extend low 8 bytes to 8 dwords.
    return _mm256_cvtepu8_epi32( m8 );
}

// Same idea, but force the immediate into the instruction stream via inline asm
// so the constant can never be hoisted to a .rodata load.
__attribute((noinline))
__m256i get_indices11()
{
    __m128i m8;
    __m256i r;
    asm("movabs $0x0706050403020100, %%rax\n\t"
        "vmovq %%rax, %0"
        : "=x"(m8) : : "rax" );
    asm("vpmovzxbd %1, %0" : "=x"(r) : "x"(m8) );
    return r;
}


#ifdef __AVX512F__
// AVX-512: 16-lane identity [0..15] as dwords in a ZMM, no memory reads.
// vpmovzxbd zmm,xmm fans 16 bytes -> 16 dwords in one op; the only extra work
// vs AVX2 is assembling the 16-byte source, which needs two movabs immediates.

// Intrinsic form -- watch whether the compiler spills the constant to .rodata.
__attribute((noinline))
__m512i get_indices512_intrin()
{
    __m128i m8 = _mm_set_epi64x( 0x0f0e0d0c0b0a0908ll, 0x0706050403020100ll );
    return _mm512_cvtepu8_epi32( m8 );
}

// Forced-immediate form -- guarantees the constants stay in the instruction stream.
__attribute((noinline))
__m512i get_indices512_asm()
{
    __m128i m8;
    __m512i r;
    asm("movabs $0x0706050403020100, %%rax\n\t"
        "movabs $0x0f0e0d0c0b0a0908, %%rdx\n\t"
        "vmovq   %%rax, %0\n\t"
        "vpinsrq $1, %%rdx, %0, %0"
        : "=x"(m8) : : "rax", "rdx" );
    asm("vpmovzxbd %1, %0" : "=v"(r) : "v"(m8) );
    return r;
}

// Carry only ONE 8-byte immediate (bytes 0..7); derive the high half via a
// broadcast + masked add in the dword domain. vpmovzxbd of the broadcasted
// qword yields [0..7, 0..7]; a merge-masked +8 fixes up the top 8 lanes.
__attribute((noinline))
__m512i get_indices512_bcastadd()
{
    __m128i m8;
    asm("movabs $0x0706050403020100, %%rax\n\t"
        "vmovq        %%rax, %0\n\t"
        "vpbroadcastq %0, %0"
        : "=x"(m8) : : "rax" );
    __m512i base = _mm512_cvtepu8_epi32( m8 );   // [0..7, 0..7]
    __m512i eights;
    asm("vpbroadcastd %1, %0" : "=v"(eights) : "r"(8) );  // all 8s, from GPR (no mem)
    return _mm512_mask_add_epi32( base, 0xFF00, base, eights ); // +8 to high lanes
}

// Derive the high 8 source bytes with vpaddb (0..7 + 8 = 8..15, no carry),
// stitch the 16 bytes with vpunpcklqdq, then a single vpmovzxbd.
__attribute((noinline))
__m512i get_indices512_addb()
{
    __m128i m8;
    asm("movabs $0x0706050403020100, %%rax\n\t"
        "vmovq %%rax, %0"
        : "=x"(m8) : : "rax" );
    __m128i eight_b;
    asm("vpbroadcastb %1, %0" : "=x"(eight_b) : "r"(8) ); // low bytes = 0x08
    __m128i hi = _mm_add_epi8( m8, eight_b );             // bytes 8..15 in low qword
    __m128i bytes16 = _mm_unpacklo_epi64( m8, hi );       // bytes 0..15
    return _mm512_cvtepu8_epi32( bytes16 );
}

// Build 0..7 once, derive 8..15 with a dword add, glue halves with vinserti32x8.
__attribute((noinline))
__m512i get_indices512_halves()
{
    __m128i m8;
    asm("movabs $0x0706050403020100, %%rax\n\t"
        "vmovq %%rax, %0"
        : "=x"(m8) : : "rax" );
    __m256i lo = _mm256_cvtepu8_epi32( m8 );   // 0..7
    __m256i eights;
    asm("vpbroadcastd %1, %0" : "=v"(eights) : "r"(8) );
    __m256i hi = _mm256_add_epi32( lo, eights );  // 8..15
    return _mm512_inserti32x8( _mm512_castsi256_si512( lo ), hi, 1 );
}

typedef __m512i (*pfnGetIndices512)(void);

double
timeFn512( pfnGetIndices512 pfn, size_t N )
{
    for ( size_t i = 0; i < N; i++ ) {
        volatile __m512i m = pfn();
    }
    uint64_t start = __rdtsc();
    for ( size_t i = 0; i < N; i++ ) {
        volatile __m512i m = pfn();
    }
    uint64_t et = __rdtsc() - start;
    return (double) et / N;
}
#endif // __AVX512F__

double
timeFn( pfnGetIndices pfn, size_t N )
{
    for ( size_t i = 0; i < N; i++ ) {
        volatile __m256i m = pfn();
    }
    uint64_t start = __rdtsc();
    for ( size_t i = 0; i < N; i++ ) {
        volatile __m256i m = pfn();
    }
    uint64_t et = __rdtsc() - start;
    return (double) et / N;
}

int
main()
{
    size_t N = 100000000;

    // Runtime designed for JIT - it holds relocated functions and controls their lifetime.
    JitRuntime rt;
    Func get_indices8;
    if ( 0 != emitIdentityIndices( rt, &get_indices8 ) ) {
        fprintf( stderr, "Error emitting custom identity index function\n" );
        exit( 1 );
    }

    printf( "GetIndices8: %.2f clocks/iteration\n", timeFn( get_indices8, N ) );

    printf( "GetIndices1: %.2f clocks/iteration\n", timeFn( get_indices1, N ) );
    printf( "GetIndices2: %.2f clocks/iteration\n", timeFn( get_indices2, N ) );
    printf( "GetIndices3: %.2f clocks/iteration\n", timeFn( get_indices3, N ) );
    printf( "GetIndices4: %.2f clocks/iteration\n", timeFn( get_indices4, N ) );
    printf( "GetIndices5: %.2f clocks/iteration\n", timeFn( get_indices5, N ) );
    printf( "GetIndices6: %.2f clocks/iteration\n", timeFn( get_indices6, N ) );
    printf( "GetIndices7: %.2f clocks/iteration\n", timeFn( get_indices7, N ) );
    printf( "GetIndices10: %.2f clocks/iteration\n", timeFn( get_indices10, N ) );
    printf( "GetIndices11: %.2f clocks/iteration\n", timeFn( get_indices11, N ) );

    print( get_indices10() );
    print( get_indices11() );
    print( get_indices1() );

#ifdef __AVX512F__
    printf( "--- AVX-512 (16-lane) ---\n" );
    printf( "GetIndices512_intrin:   %.2f clocks/iteration\n", timeFn512( get_indices512_intrin, N ) );
    printf( "GetIndices512_asm:      %.2f clocks/iteration\n", timeFn512( get_indices512_asm, N ) );
    printf( "GetIndices512_bcastadd: %.2f clocks/iteration\n", timeFn512( get_indices512_bcastadd, N ) );
    printf( "GetIndices512_addb:     %.2f clocks/iteration\n", timeFn512( get_indices512_addb, N ) );
    printf( "GetIndices512_halves:   %.2f clocks/iteration\n", timeFn512( get_indices512_halves, N ) );
    print( get_indices512_asm() );
    print( get_indices512_bcastadd() );
    print( get_indices512_addb() );
    print( get_indices512_halves() );
#endif
    print( get_indices2() );
    print( get_indices3() );
    print( get_indices4() );
    print( get_indices5() );
    print( get_indices8() );
    return 0;
}

