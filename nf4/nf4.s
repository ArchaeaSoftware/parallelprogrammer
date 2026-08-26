	.text
	.file	"nf4.cpp"
	.section	.rodata.cst4,"aM",@progbits,4
	.p2align	2               # -- Begin function _Z14dQuantizeNF4_0f
.LCPI0_0:
	.long	1025702655              # float 0.0397901498
.LCPI0_1:
	.long	3199068790              # float -0.33967942
.LCPI0_2:
	.long	3206304368              # float -0.610632896
.LCPI0_3:
	.long	3210288345              # float -0.84809643
.LCPI0_4:
	.long	3203105920              # float -0.45999527
.LCPI0_5:
	.long	3188537532              # float -0.137911737
.LCPI0_6:
	.long	3195026668              # float -0.234607399
.LCPI0_7:
	.long	3174725745              # float -0.0455250181
.LCPI0_8:
	.long	1053250553              # float 0.389312536
.LCPI0_9:
	.long	1045456864              # float 0.203521252
.LCPI0_10:
	.long	1039550563              # float 0.120255254
.LCPI0_11:
	.long	1049985748              # float 0.292013764
.LCPI0_12:
	.long	1059360175              # float 0.64278692
.LCPI0_13:
	.long	1056992515              # float 0.501663387
.LCPI0_14:
	.long	1063029209              # float 0.861478388
	.text
	.globl	_Z14dQuantizeNF4_0f
	.p2align	4, 0x90
	.type	_Z14dQuantizeNF4_0f,@function
_Z14dQuantizeNF4_0f:                    # @_Z14dQuantizeNF4_0f
	.cfi_startproc
# %bb.0:
	vucomiss	.LCPI0_0(%rip), %xmm0
	jbe	.LBB0_8
# %bb.1:
	vucomiss	.LCPI0_8(%rip), %xmm0
	jbe	.LBB0_5
# %bb.2:
	vucomiss	.LCPI0_12(%rip), %xmm0
	jbe	.LBB0_4
# %bb.3:
	vucomiss	.LCPI0_14(%rip), %xmm0
	seta	%al
	orb	$14, %al
	retq
.LBB0_8:
	vucomiss	.LCPI0_1(%rip), %xmm0
	jbe	.LBB0_12
# %bb.9:
	vucomiss	.LCPI0_5(%rip), %xmm0
	jbe	.LBB0_11
# %bb.10:
	vucomiss	.LCPI0_7(%rip), %xmm0
	seta	%al
	orb	$6, %al
	retq
.LBB0_5:
	vucomiss	.LCPI0_9(%rip), %xmm0
	jbe	.LBB0_7
# %bb.6:
	vucomiss	.LCPI0_11(%rip), %xmm0
	seta	%al
	orb	$10, %al
	retq
.LBB0_12:
	vucomiss	.LCPI0_2(%rip), %xmm0
	jbe	.LBB0_14
# %bb.13:
	vucomiss	.LCPI0_4(%rip), %xmm0
	seta	%al
	orb	$2, %al
	retq
.LBB0_4:
	vucomiss	.LCPI0_13(%rip), %xmm0
	seta	%al
	orb	$12, %al
	retq
.LBB0_11:
	vucomiss	.LCPI0_6(%rip), %xmm0
	seta	%al
	orb	$4, %al
	retq
.LBB0_7:
	vucomiss	.LCPI0_10(%rip), %xmm0
	seta	%al
	orb	$8, %al
	retq
.LBB0_14:
	vucomiss	.LCPI0_3(%rip), %xmm0
	seta	%al
	retq
.Lfunc_end0:
	.size	_Z14dQuantizeNF4_0f, .Lfunc_end0-_Z14dQuantizeNF4_0f
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst4,"aM",@progbits,4
	.p2align	2               # -- Begin function _Z14float_to_NF4_0PKfmf
.LCPI1_0:
	.long	2147483647              # float NaN
	.text
	.globl	_Z14float_to_NF4_0PKfmf
	.p2align	4, 0x90
	.type	_Z14float_to_NF4_0PKfmf,@function
_Z14float_to_NF4_0PKfmf:                # @_Z14float_to_NF4_0PKfmf
	.cfi_startproc
# %bb.0:
	vmovsd	(%rdi), %xmm1           # xmm1 = mem[0],zero
	vbroadcastss	%xmm0, %xmm2
	vsubps	%xmm2, %xmm1, %xmm2
	vbroadcastss	.LCPI1_0(%rip), %xmm1 # xmm1 = [NaN,NaN,NaN,NaN]
	vandps	%xmm1, %xmm2, %xmm2
	vpermilps	$245, %xmm2, %xmm3 # xmm3 = xmm2[1,1,3,3]
	xorl	%eax, %eax
	vucomiss	%xmm3, %xmm2
	seta	%al
	vminss	%xmm2, %xmm3, %xmm2
	vmovss	8(%rdi), %xmm3          # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$2, %ecx
	cmovbel	%eax, %ecx
	vmovss	12(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	movl	$3, %eax
	cmovbel	%ecx, %eax
	vminss	%xmm2, %xmm3, %xmm2
	vmovss	16(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$4, %ecx
	cmovbel	%eax, %ecx
	vmovss	20(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$5, %eax
	cmovbel	%ecx, %eax
	vmovss	24(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	movl	$6, %ecx
	cmovbel	%eax, %ecx
	vminss	%xmm2, %xmm3, %xmm2
	vmovss	28(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$7, %eax
	cmovbel	%ecx, %eax
	vmovss	32(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$8, %ecx
	cmovbel	%eax, %ecx
	vmovss	36(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	movl	$9, %eax
	cmovbel	%ecx, %eax
	vminss	%xmm2, %xmm3, %xmm2
	vmovss	40(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$10, %ecx
	cmovbel	%eax, %ecx
	vmovss	44(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$11, %eax
	cmovbel	%ecx, %eax
	vmovss	48(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	movl	$12, %ecx
	cmovbel	%eax, %ecx
	vminss	%xmm2, %xmm3, %xmm2
	vmovss	52(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$13, %eax
	cmovbel	%ecx, %eax
	vmovss	56(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm3
	vandps	%xmm1, %xmm3, %xmm3
	vucomiss	%xmm3, %xmm2
	vminss	%xmm2, %xmm3, %xmm2
	movl	$14, %ecx
	cmovbel	%eax, %ecx
	vmovss	60(%rdi), %xmm3         # xmm3 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm3, %xmm0
	vandps	%xmm1, %xmm0, %xmm0
	vucomiss	%xmm0, %xmm2
	movl	$15, %eax
	cmovbel	%ecx, %eax
	retq
.Lfunc_end1:
	.size	_Z14float_to_NF4_0PKfmf, .Lfunc_end1-_Z14float_to_NF4_0PKfmf
	.cfi_endproc
                                        # -- End function
	.globl	_Z16float_to_NF4_midPKfmf # -- Begin function _Z16float_to_NF4_midPKfmf
	.p2align	4, 0x90
	.type	_Z16float_to_NF4_midPKfmf,@function
_Z16float_to_NF4_midPKfmf:              # @_Z16float_to_NF4_midPKfmf
	.cfi_startproc
# %bb.0:
	xorl	%eax, %eax
	vucomiss	32(%rdi), %xmm0
	setae	%al
	leal	4(,%rax,8), %ecx
                                        # kill: def $eax killed $eax killed $rax def $rax
	shll	$3, %eax
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	leal	2(%rax), %ecx
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	leal	1(%rax), %ecx
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
                                        # kill: def $eax killed $eax killed $rax
	retq
.Lfunc_end2:
	.size	_Z16float_to_NF4_midPKfmf, .Lfunc_end2-_Z16float_to_NF4_midPKfmf
	.cfi_endproc
                                        # -- End function
	.globl	_Z14float_to_NF4_1PKfmf # -- Begin function _Z14float_to_NF4_1PKfmf
	.p2align	4, 0x90
	.type	_Z14float_to_NF4_1PKfmf,@function
_Z14float_to_NF4_1PKfmf:                # @_Z14float_to_NF4_1PKfmf
	.cfi_startproc
# %bb.0:
	xorl	%eax, %eax
	vucomiss	32(%rdi), %xmm0
	setae	%al
	leal	4(,%rax,8), %ecx
                                        # kill: def $eax killed $eax killed $rax def $rax
	shll	$3, %eax
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	leal	2(%rax), %ecx
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	leal	1(%rax), %ecx
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	cltq
	vmovss	4(%rdi,%rax,4), %xmm1   # xmm1 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm1, %xmm1
	vsubss	(%rdi,%rax,4), %xmm0, %xmm0
	xorl	%ecx, %ecx
	vucomiss	%xmm1, %xmm0
	seta	%cl
	addl	%ecx, %eax
                                        # kill: def $eax killed $eax killed $rax
	retq
.Lfunc_end3:
	.size	_Z14float_to_NF4_1PKfmf, .Lfunc_end3-_Z14float_to_NF4_1PKfmf
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst4,"aM",@progbits,4
	.p2align	2               # -- Begin function _Z15float_to_NF4_16PjPKfmS1_
.LCPI4_0:
	.long	8                       # 0x8
.LCPI4_1:
	.long	4                       # 0x4
.LCPI4_2:
	.long	2                       # 0x2
.LCPI4_3:
	.long	1                       # 0x1
	.text
	.globl	_Z15float_to_NF4_16PjPKfmS1_
	.p2align	4, 0x90
	.type	_Z15float_to_NF4_16PjPKfmS1_,@function
_Z15float_to_NF4_16PjPKfmS1_:           # @_Z15float_to_NF4_16PjPKfmS1_
	.cfi_startproc
# %bb.0:
	vbroadcastss	_ZL11NF4_LUT_mid+32(%rip), %zmm0
	vmovaps	(%rcx), %zmm1
	vcmpleps	%zmm1, %zmm0, %k1
	vpbroadcastd	.LCPI4_0(%rip), %zmm0 {%k1} {z}
	vmovdqa64	_ZL11NF4_LUT_mid(%rip), %zmm2
	vpbroadcastd	.LCPI4_1(%rip), %zmm3 # zmm3 = [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]
	vpord	%zmm3, %zmm0, %zmm4
	vpermd	%zmm2, %zmm4, %zmm4
	vcmpleps	%zmm1, %zmm4, %k1
	vpbroadcastd	.LCPI4_2(%rip), %zmm4 # zmm4 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
	vpord	%zmm3, %zmm0, %zmm0 {%k1}
	vpord	%zmm4, %zmm0, %zmm3
	vpermd	%zmm2, %zmm3, %zmm3
	vcmpleps	%zmm1, %zmm3, %k1
	vpord	%zmm4, %zmm0, %zmm0 {%k1}
	vpbroadcastd	.LCPI4_3(%rip), %zmm3 # zmm3 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
	vpord	%zmm3, %zmm0, %zmm4
	vpermd	%zmm2, %zmm4, %zmm2
	vcmpleps	%zmm1, %zmm2, %k1
	vpord	%zmm3, %zmm0, %zmm0 {%k1}
	vmovdqa64	%zmm0, (%rdi)
	vzeroupper
	retq
.Lfunc_end4:
	.size	_Z15float_to_NF4_16PjPKfmS1_, .Lfunc_end4-_Z15float_to_NF4_16PjPKfmS1_
	.cfi_endproc
                                        # -- End function
	.globl	_Z14float_to_NF4_2PKfmf # -- Begin function _Z14float_to_NF4_2PKfmf
	.p2align	4, 0x90
	.type	_Z14float_to_NF4_2PKfmf,@function
_Z14float_to_NF4_2PKfmf:                # @_Z14float_to_NF4_2PKfmf
	.cfi_startproc
# %bb.0:
	vxorps	%xmm1, %xmm1, %xmm1
	xorl	%eax, %eax
	vucomiss	%xmm0, %xmm1
	jbe	.LBB5_2
# %bb.1:
	vucomiss	16(%rdi), %xmm0
	setae	%al
	leal	2(,%rax,4), %ecx
                                        # kill: def $eax killed $eax killed $rax def $rax
	shll	$2, %eax
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	leal	1(%rax), %ecx
	vucomiss	(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	jmp	.LBB5_3
.LBB5_2:
	vucomiss	44(%rdi), %xmm0
	setae	%al
	leal	2(,%rax,4), %ecx
                                        # kill: def $eax killed $eax killed $rax def $rax
	shll	$2, %eax
	vucomiss	28(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	leal	1(%rax), %ecx
	vucomiss	28(%rdi,%rcx,4), %xmm0
	cmovael	%ecx, %eax
	addl	$7, %eax
.LBB5_3:
	cltq
	vmovss	4(%rdi,%rax,4), %xmm1   # xmm1 = mem[0],zero,zero,zero
	vsubss	%xmm0, %xmm1, %xmm1
	vsubss	(%rdi,%rax,4), %xmm0, %xmm0
	xorl	%ecx, %ecx
	vucomiss	%xmm1, %xmm0
	seta	%cl
	addl	%ecx, %eax
                                        # kill: def $eax killed $eax killed $rax
	retq
.Lfunc_end5:
	.size	_Z14float_to_NF4_2PKfmf, .Lfunc_end5-_Z14float_to_NF4_2PKfmf
	.cfi_endproc
                                        # -- End function
	.globl	_Z14float_to_NF4_3PKfmf # -- Begin function _Z14float_to_NF4_3PKfmf
	.p2align	4, 0x90
	.type	_Z14float_to_NF4_3PKfmf,@function
_Z14float_to_NF4_3PKfmf:                # @_Z14float_to_NF4_3PKfmf
	.cfi_startproc
# %bb.0:
	pushq	%rax
	.cfi_def_cfa_offset 16
	callq	_Z14dQuantizeNF4_0f
	movzbl	%al, %eax
	popq	%rcx
	.cfi_def_cfa_offset 8
	retq
.Lfunc_end6:
	.size	_Z14float_to_NF4_3PKfmf, .Lfunc_end6-_Z14float_to_NF4_3PKfmf
	.cfi_endproc
                                        # -- End function
	.section	.rodata.cst4,"aM",@progbits,4
	.p2align	2               # -- Begin function main
.LCPI7_0:
	.long	805306368               # float 4.65661287E-10
.LCPI7_1:
	.long	1065353216              # float 1
.LCPI7_2:
	.long	2147483648              # float -0
.LCPI7_4:
	.long	1034091263              # float 0.0795802995
.LCPI7_5:
	.long	8                       # 0x8
.LCPI7_6:
	.long	4                       # 0x4
.LCPI7_7:
	.long	2                       # 0x2
.LCPI7_8:
	.long	1                       # 0x1
	.section	.rodata.cst8,"aM",@progbits,8
	.p2align	3
.LCPI7_3:
	.quad	4485585228861014016     # double 7.4505805969238281E-9
	.text
	.globl	main
	.p2align	4, 0x90
	.type	main,@function
main:                                   # @main
	.cfi_startproc
# %bb.0:
	pushq	%rbp
	.cfi_def_cfa_offset 16
	pushq	%rbx
	.cfi_def_cfa_offset 24
	subq	$232, %rsp
	.cfi_def_cfa_offset 256
	.cfi_offset %rbx, -24
	.cfi_offset %rbp, -16
	xorl	%ebx, %ebx
	vbroadcastss	.LCPI7_2(%rip), %xmm0 # xmm0 = [-0.0E+0,-0.0E+0,-0.0E+0,-0.0E+0]
	vmovaps	%xmm0, 64(%rsp)         # 16-byte Spill
	.p2align	4, 0x90
.LBB7_1:                                # =>This Inner Loop Header: Depth=1
	callq	rand
	vcvtsi2ss	%eax, %xmm2, %xmm0
	vmulss	.LCPI7_0(%rip), %xmm0, %xmm0
	vucomiss	.LCPI7_1(%rip), %xmm0
	jne	.LBB7_2
	jnp	.LBB7_1
.LBB7_2:                                #   in Loop: Header=BB7_1 Depth=1
	vmovaps	%xmm0, (%rsp)           # 16-byte Spill
	callq	rand
	vmovaps	(%rsp), %xmm1           # 16-byte Reload
	vxorps	64(%rsp), %xmm1, %xmm0  # 16-byte Folded Reload
	kmovw	%eax, %k1
	vmovss	%xmm0, %xmm1, %xmm1 {%k1}
	vmovss	%xmm1, _ZZ4mainE2in(,%rbx,4)
	addq	$1, %rbx
	cmpq	$134217728, %rbx        # imm = 0x8000000
	jne	.LBB7_1
# %bb.3:
	rdtsc
	movq	%rdx, %rbx
	shlq	$32, %rbx
	orq	%rax, %rbx
	xorl	%ebp, %ebp
	.p2align	4, 0x90
.LBB7_4:                                # =>This Inner Loop Header: Depth=1
	vmovss	_ZZ4mainE2in(%rbp), %xmm0 # xmm0 = mem[0],zero,zero,zero
	movl	$_ZL7NF4_LUT, %edi
	callq	_Z14float_to_NF4_0PKfmf
	movl	%eax, _ZZ4mainE4inx0(%rbp)
	addq	$4, %rbp
	cmpq	$536870912, %rbp        # imm = 0x20000000
	jne	.LBB7_4
# %bb.5:
	rdtsc
	shlq	$32, %rdx
	orq	%rax, %rdx
	subq	%rbx, %rdx
	vcvtusi2sd	%rdx, %xmm2, %xmm0
	vmulsd	.LCPI7_3(%rip), %xmm0, %xmm0
	movl	$.L.str, %edi
	movb	$1, %al
	callq	printf
	rdtsc
	movq	%rdx, %rcx
	shlq	$32, %rcx
	orq	%rax, %rcx
	vbroadcastss	.LCPI7_4(%rip), %zmm0 # zmm0 = [7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2,7.95802995E-2]
	vpbroadcastd	.LCPI7_6(%rip), %zmm8 # zmm8 = [4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]
	vpbroadcastd	.LCPI7_7(%rip), %zmm9 # zmm9 = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2]
	xorl	%eax, %eax
	vpbroadcastd	.LCPI7_8(%rip), %zmm10 # zmm10 = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
	vpternlogd	$255, %zmm1, %zmm1, %zmm1
	.p2align	4, 0x90
.LBB7_6:                                # =>This Inner Loop Header: Depth=1
	vmovaps	_ZZ4mainE2in(%rax), %zmm2
	vcmpgeps	%zmm0, %zmm2, %k1
	vpbroadcastd	.LCPI7_5(%rip), %zmm3 {%k1} {z}
	vpord	%zmm8, %zmm3, %zmm4
	kxnorw	%k0, %k0, %k1
	vgatherdps	_ZL7NF4_LUT(,%zmm4,4), %zmm5 {%k1}
	vcmpnleps	%zmm2, %zmm5, %k1
	vmovdqa32	%zmm3, %zmm4 {%k1}
	vpord	%zmm9, %zmm4, %zmm3
	kxnorw	%k0, %k0, %k1
	vgatherdps	_ZL7NF4_LUT(,%zmm3,4), %zmm5 {%k1}
	vcmpnleps	%zmm2, %zmm5, %k1
	vmovdqa32	%zmm4, %zmm3 {%k1}
	vpord	%zmm10, %zmm3, %zmm4
	vpmovzxdq	%ymm4, %zmm5    # zmm5 = ymm4[0],zero,ymm4[1],zero,ymm4[2],zero,ymm4[3],zero,ymm4[4],zero,ymm4[5],zero,ymm4[6],zero,ymm4[7],zero
	vextracti64x4	$1, %zmm4, %ymm6
	vpmovzxdq	%ymm6, %zmm6    # zmm6 = ymm6[0],zero,ymm6[1],zero,ymm6[2],zero,ymm6[3],zero,ymm6[4],zero,ymm6[5],zero,ymm6[6],zero,ymm6[7],zero
	kxnorw	%k0, %k0, %k1
	vgatherqps	_ZL7NF4_LUT(,%zmm6,4), %ymm7 {%k1}
	kxnorw	%k0, %k0, %k1
	vgatherqps	_ZL7NF4_LUT(,%zmm5,4), %ymm6 {%k1}
	vinsertf64x4	$1, %ymm7, %zmm6, %zmm5
	vcmpnleps	%zmm2, %zmm5, %k1
	vmovdqa32	%zmm3, %zmm4 {%k1}
	vpsubd	%zmm1, %zmm4, %zmm3
	kxnorw	%k0, %k0, %k1
	vgatherdps	_ZL7NF4_LUT(,%zmm3,4), %zmm5 {%k1}
	vsubps	%zmm2, %zmm5, %zmm3
	kxnorw	%k0, %k0, %k1
	vgatherdps	_ZL7NF4_LUT(,%zmm4,4), %zmm5 {%k1}
	vsubps	%zmm5, %zmm2, %zmm2
	vcmpltps	%zmm2, %zmm3, %k0
	vpmovm2d	%k0, %zmm2
	vpsubd	%zmm2, %zmm4, %zmm2
	vmovdqa64	%zmm2, _ZZ4mainE4inx0(%rax)
	addq	$64, %rax
	cmpq	$536870912, %rax        # imm = 0x20000000
	jne	.LBB7_6
# %bb.7:
	vmovdqu64	%zmm10, 128(%rsp) # 64-byte Spill
	vmovdqu64	%zmm9, 64(%rsp) # 64-byte Spill
	vmovdqu64	%zmm8, (%rsp)   # 64-byte Spill
	rdtsc
	shlq	$32, %rdx
	orq	%rax, %rdx
	subq	%rcx, %rdx
	vcvtusi2sd	%rdx, %xmm11, %xmm0
	vmulsd	.LCPI7_3(%rip), %xmm0, %xmm0
	movl	$.L.str.1, %edi
	movb	$1, %al
	vzeroupper
	callq	printf
	rdtsc
	movq	%rdx, %rbx
	shlq	$32, %rbx
	orq	%rax, %rbx
	xorl	%ebp, %ebp
	.p2align	4, 0x90
.LBB7_8:                                # =>This Inner Loop Header: Depth=1
	vmovss	_ZZ4mainE2in(%rbp), %xmm0 # xmm0 = mem[0],zero,zero,zero
	callq	_Z14dQuantizeNF4_0f
	movzbl	%al, %eax
	movl	%eax, _ZZ4mainE4inx1(%rbp)
	addq	$4, %rbp
	cmpq	$536870912, %rbp        # imm = 0x20000000
	jne	.LBB7_8
# %bb.9:
	rdtsc
	shlq	$32, %rdx
	orq	%rax, %rdx
	subq	%rbx, %rdx
	vcvtusi2sd	%rdx, %xmm11, %xmm0
	vmulsd	.LCPI7_3(%rip), %xmm0, %xmm0
	movl	$.L.str.2, %edi
	movb	$1, %al
	callq	printf
	rdtsc
	movq	%rdx, %rcx
	shlq	$32, %rcx
	orq	%rax, %rcx
	vmovdqa64	_ZL11NF4_LUT_mid(%rip), %zmm0
	vbroadcastss	_ZL11NF4_LUT_mid+32(%rip), %zmm1
	xorl	%eax, %eax
	vmovdqu64	(%rsp), %zmm5   # 64-byte Reload
	vmovdqu64	64(%rsp), %zmm6 # 64-byte Reload
	vmovdqu64	128(%rsp), %zmm7 # 64-byte Reload
	.p2align	4, 0x90
.LBB7_10:                               # =>This Inner Loop Header: Depth=1
	vmovaps	_ZZ4mainE2in(%rax), %zmm2
	vcmpleps	%zmm2, %zmm1, %k1
	vpbroadcastd	.LCPI7_5(%rip), %zmm3 {%k1} {z}
	vpord	%zmm5, %zmm3, %zmm4
	vpermd	%zmm0, %zmm4, %zmm4
	vcmpleps	%zmm2, %zmm4, %k1
	vpord	%zmm5, %zmm3, %zmm3 {%k1}
	vpord	%zmm6, %zmm3, %zmm4
	vpermd	%zmm0, %zmm4, %zmm4
	vcmpleps	%zmm2, %zmm4, %k1
	vpord	%zmm6, %zmm3, %zmm3 {%k1}
	vpord	%zmm7, %zmm3, %zmm4
	vpermd	%zmm0, %zmm4, %zmm4
	vcmpleps	%zmm2, %zmm4, %k1
	vpord	%zmm7, %zmm3, %zmm3 {%k1}
	vmovdqa64	%zmm3, _ZZ4mainE4inx1(%rax)
	addq	$64, %rax
	cmpq	$536870912, %rax        # imm = 0x20000000
	jne	.LBB7_10
# %bb.11:
	rdtsc
	shlq	$32, %rdx
	orq	%rax, %rdx
	subq	%rcx, %rdx
	vcvtusi2sd	%rdx, %xmm11, %xmm0
	vmulsd	.LCPI7_3(%rip), %xmm0, %xmm0
	movl	$.L.str.3, %edi
	movb	$1, %al
	vzeroupper
	callq	printf
	xorl	%ebx, %ebx
	xorl	%ebp, %ebp
	jmp	.LBB7_12
	.p2align	4, 0x90
.LBB7_14:                               #   in Loop: Header=BB7_12 Depth=1
	addq	$4, %rbx
	cmpq	$536870912, %rbx        # imm = 0x20000000
	je	.LBB7_15
.LBB7_12:                               # =>This Inner Loop Header: Depth=1
	movl	_ZZ4mainE4inx0(%rbx), %edx
	movl	_ZZ4mainE4inx1(%rbx), %esi
	cmpl	%esi, %edx
	je	.LBB7_14
# %bb.13:                               #   in Loop: Header=BB7_12 Depth=1
	movl	$.L.str.4, %edi
	xorl	%eax, %eax
	callq	printf
	vmovss	_ZZ4mainE2in(%rbx), %xmm0 # xmm0 = mem[0],zero,zero,zero
	vmovss	%xmm0, (%rsp)           # 4-byte Spill
	movl	$_ZL7NF4_LUT, %edi
	callq	_Z14float_to_NF4_0PKfmf
	movl	%eax, _ZZ4mainE4inx0(%rbx)
	vmovss	(%rsp), %xmm0           # 4-byte Reload
                                        # xmm0 = mem[0],zero,zero,zero
	callq	_Z14dQuantizeNF4_0f
	movzbl	%al, %eax
	movl	%eax, _ZZ4mainE4inx1(%rbx)
	addl	$1, %ebp
	jmp	.LBB7_14
.LBB7_15:
	movl	$.L.str.5, %edi
	movl	%ebp, %esi
	xorl	%eax, %eax
	callq	printf
	xorl	%eax, %eax
	addq	$232, %rsp
	.cfi_def_cfa_offset 24
	popq	%rbx
	.cfi_def_cfa_offset 16
	popq	%rbp
	.cfi_def_cfa_offset 8
	retq
.Lfunc_end7:
	.size	main, .Lfunc_end7-main
	.cfi_endproc
                                        # -- End function
	.type	_ZL11NF4_LUT_mid,@object # @_ZL11NF4_LUT_mid
	.section	.rodata,"a",@progbits
	.p2align	6
_ZL11NF4_LUT_mid:
	.long	0                       # float 0
	.long	3210288344              # float -0.84809637
	.long	3206304368              # float -0.610632896
	.long	3203105920              # float -0.45999527
	.long	3199068790              # float -0.33967942
	.long	3195026668              # float -0.234607399
	.long	3188537532              # float -0.137911737
	.long	3174725745              # float -0.0455250181
	.long	1025702655              # float 0.0397901498
	.long	1039550562              # float 0.120255247
	.long	1045456864              # float 0.203521252
	.long	1049985748              # float 0.292013764
	.long	1053250553              # float 0.389312536
	.long	1056992516              # float 0.501663446
	.long	1059360175              # float 0.64278692
	.long	1063029210              # float 0.861478447
	.size	_ZL11NF4_LUT_mid, 64

	.type	_ZZ4mainE2in,@object    # @_ZZ4mainE2in
	.local	_ZZ4mainE2in
	.comm	_ZZ4mainE2in,536870912,64
	.type	_ZZ4mainE4inx0,@object  # @_ZZ4mainE4inx0
	.local	_ZZ4mainE4inx0
	.comm	_ZZ4mainE4inx0,536870912,64
	.type	_ZZ4mainE4inx1,@object  # @_ZZ4mainE4inx1
	.local	_ZZ4mainE4inx1
	.comm	_ZZ4mainE4inx1,536870912,64
	.type	_ZL7NF4_LUT,@object     # @_ZL7NF4_LUT
	.p2align	4
_ZL7NF4_LUT:
	.long	3212836864              # float -1
	.long	3207739825              # float -0.6961928
	.long	3204868912              # float -0.525073051
	.long	3200922272              # float -0.394917488
	.long	3197215309              # float -0.284441382
	.long	3191682367              # float -0.18477343
	.long	3183114353              # float -0.0910500362
	.long	0                       # float 0
	.long	1034091263              # float 0.0795802995
	.long	1042598627              # float 0.160930201
	.long	1048315101              # float 0.246112302
	.long	1051525946              # float 0.337915242
	.long	1054975160              # float 0.440709829
	.long	1058015147              # float 0.562617004
	.long	1060705203              # float 0.722956836
	.long	1065353216              # float 1
	.size	_ZL7NF4_LUT, 64

	.type	.L.str,@object          # @.str
	.section	.rodata.str1.1,"aMS",@progbits,1
.L.str:
	.asciz	"gold: %.2f clocks/iteration\n"
	.size	.L.str, 29

	.type	.L.str.1,@object        # @.str.1
.L.str.1:
	.asciz	"binsearch: %.2f clocks/iteration\n"
	.size	.L.str.1, 34

	.type	.L.str.2,@object        # @.str.2
.L.str.2:
	.asciz	"bitsandbytes: %.2f clocks/iteration\n"
	.size	.L.str.2, 37

	.type	.L.str.3,@object        # @.str.3
.L.str.3:
	.asciz	"AVX512: %.2f clocks/iteration\n"
	.size	.L.str.3, 31

	.type	.L.str.4,@object        # @.str.4
.L.str.4:
	.asciz	"Mismatch: %u should be %u\n"
	.size	.L.str.4, 27

	.type	.L.str.5,@object        # @.str.5
.L.str.5:
	.asciz	"%d errors\n"
	.size	.L.str.5, 11

	.ident	"clang version 10.0.0-4ubuntu1 "
	.section	".note.GNU-stack","",@progbits
	.addrsig
	.addrsig_sym _ZL7NF4_LUT
