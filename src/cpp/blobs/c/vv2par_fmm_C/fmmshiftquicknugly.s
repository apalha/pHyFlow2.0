	.file	"fmmshiftquicknugly.cpp"
	.section	.ctors,"w"
	.align 4
	.long	__GLOBAL__I__Z11shift_allociii
.lcomm _I,16
.lcomm _rminshift,16
.lcomm _pshift,16
.lcomm _wksp,16
.lcomm _potshift,16
.lcomm _rpow,16
	.section .rdata,"dr"
	.align 4
LC2:
	.long	-998637568
	.text
	.align 2
	.p2align 4,,15
.globl __Z11shift_allociii
	.def	__Z11shift_allociii;	.scl	2;	.type	32;	.endef
__Z11shift_allociii:
	pushl	%ebp
	movl	%esp, %ebp
	pushl	%esi
	pushl	%ebx
	subl	$16, %esp
	movl	8(%ebp), %ebx
	cmpl	%ebx, _pshift
	movl	12(%ebp), %esi
	je	L1
	movl	_wksp, %edx
	movl	%edx, (%esp)
	call	_mxFree
	movl	%ebx, %eax
	sall	$5, %eax
	addl	$32, %eax
	movl	%eax, (%esp)
	call	_mxMalloc
	movl	%eax, (%esp)
	movl	%eax, _wksp
	call	_mexMakeMemoryPersistent
	movl	_rpow, %ecx
	movl	%ecx, (%esp)
	call	_mxFree
	movl	%ebx, %edx
	sall	$4, %edx
	addl	$16, %edx
	movl	%edx, (%esp)
	call	_mxMalloc
	movl	%eax, _rpow
	fld1
	fstpl	(%eax)
	fldz
	fstpl	8(%eax)
	movl	%eax, (%esp)
	call	_mexMakeMemoryPersistent
	flds	LC2
	movl	%ebx, _pshift
	leal	1(%ebx), %eax
	pushl	%eax
	fildl	(%esp)
	addl	$4, %esp
	movl	%esi, _potshift
	fdivrp	%st, %st(1)
	fstpl	(%esp)
	call	_exp2
	fstpl	_rminshift
	movl	$__Z12shift_atExitv, 8(%ebp)
	addl	$16, %esp
	popl	%ebx
	popl	%esi
	popl	%ebp
	jmp	_mexAtExit
	.p2align 4,,7
L1:
	addl	$16, %esp
	popl	%ebx
	popl	%esi
	popl	%ebp
	ret
	.align 2
	.p2align 4,,15
.globl __Z12shift_atExitv
	.def	__Z12shift_atExitv;	.scl	2;	.type	32;	.endef
__Z12shift_atExitv:
	pushl	%ebp
	movl	%esp, %ebp
	subl	$8, %esp
	movl	_rpow, %ecx
	movl	%ecx, (%esp)
	call	_mxFree
	movl	_wksp, %edx
	movl	%edx, (%esp)
	call	_mxFree
	xorl	%eax, %eax
	movl	%eax, _pshift
	leave
	ret
	.align 2
	.p2align 4,,15
.globl __Z9shift_p2pSt7complexIdES0_S0_S0_PS0_PKS0_
	.def	__Z9shift_p2pSt7complexIdES0_S0_S0_PS0_PKS0_;	.scl	2;	.type	32;	.endef
__Z9shift_p2pSt7complexIdES0_S0_S0_PS0_PKS0_:
	pushl	%ebp
	movl	$3, %edx
	movl	%esp, %ebp
	pushl	%edi
	leal	-88(%ebp), %eax
	fldz
	pushl	%esi
	pushl	%ebx
	subl	$380, %esp
	cmpl	$-1, %edx
	je	L8
L151:
	movl	$0, (%eax)
	decl	%edx
	movl	$0, 4(%eax)
	movl	$0, 8(%eax)
	movl	$0, 12(%eax)
	addl	$16, %eax
	cmpl	$-1, %edx
	jne	L151
L8:
	fldl	8(%ebp)
	xorl	%eax, %eax
	movl	$0, -104(%ebp)
	movl	_pshift, %ebx
	movl	$0, -100(%ebp)
	movl	$0, -96(%ebp)
	fstpl	-88(%ebp)
	fldl	16(%ebp)
	movl	$0, -92(%ebp)
	movl	$0, -120(%ebp)
	movl	$0, -116(%ebp)
	fstpl	-80(%ebp)
	fldl	24(%ebp)
	movl	$0, -112(%ebp)
	movl	$0, -108(%ebp)
	movl	%eax, -348(%ebp)
	fstpl	-72(%ebp)
	fldl	32(%ebp)
	fstpl	-64(%ebp)
	fldl	40(%ebp)
	fstpl	-56(%ebp)
	fldl	48(%ebp)
	fstpl	-48(%ebp)
	fldl	56(%ebp)
	fstpl	-40(%ebp)
	fldl	64(%ebp)
	fstpl	-32(%ebp)
L118:
	movl	-348(%ebp), %edx
	leal	-24(%ebp), %ecx
	movl	76(%ebp), %edi
	sall	$4, %edx
	leal	(%edx,%ecx), %esi
	movl	$1, %ecx
	fldl	-64(%esi)
	cmpl	%ebx, %ecx
	fstl	-104(%ebp)
	fldl	-56(%esi)
	movl	_wksp, %esi
	fstl	-96(%ebp)
	fldl	(%edi)
	fldl	8(%edi)
	fld	%st(1)
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(3)
	fmulp	%st, %st(4)
	fxch	%st(4)
	fmulp	%st, %st(1)
	fxch	%st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(2)
	fstl	-152(%ebp)
	fstl	-136(%ebp)
	fxch	%st(1)
	fstl	-144(%ebp)
	fstpl	-128(%ebp)
	fstpl	(%esi)
	fldl	-128(%ebp)
	fstpl	8(%esi)
	jg	L122
	leal	-88(%ebp), %eax
	addl	%edx, %eax
	movl	%eax, -364(%ebp)
	.p2align 4,,15
L31:
	fldl	-96(%ebp)
	movl	-364(%ebp), %edx
	movl	76(%ebp), %edi
	fldl	-104(%ebp)
	fld	%st(1)
	fldl	8(%edx)
	fld	%st(2)
	fldl	(%edx)
	fxch	%st(3)
	fmul	%st(2), %st
	fxch	%st(1)
	movl	%ecx, %edx
	sall	$4, %edx
	incl	%ecx
	fmul	%st(3), %st
	fxch	%st(3)
	leal	(%edx,%edi), %eax
	cmpl	%ebx, %ecx
	fmulp	%st, %st(5)
	fxch	%st(3)
	fmulp	%st, %st(1)
	fxch	%st(1)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstl	-104(%ebp)
	fxch	%st(1)
	fstl	-96(%ebp)
	fldl	(%eax)
	fldl	8(%eax)
	fld	%st(1)
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(3)
	fmulp	%st, %st(4)
	fxch	%st(4)
	fmulp	%st, %st(1)
	fxch	%st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(2)
	fstl	-184(%ebp)
	fstl	-168(%ebp)
	fxch	%st(1)
	fstl	-176(%ebp)
	fstpl	-160(%ebp)
	fstpl	(%edx,%esi)
	fldl	-160(%ebp)
	fstpl	8(%edx,%esi)
	jle	L31
L122:
	xorl	%ecx, %ecx
	cmpl	%ebx, %ecx
	jg	L124
	.p2align 4,,15
L152:
	movl	%ebx, %edx
	subl	%ecx, %edx
	cmpl	%ebx, %edx
	jge	L126
	movl	%edx, %eax
	sall	$4, %eax
	addl	%esi, %eax
	.p2align 4,,15
L41:
	fldl	16(%eax)
	incl	%edx
	fsubrl	(%eax)
	fstpl	(%eax)
	fldl	24(%eax)
	fsubrl	8(%eax)
	fstpl	8(%eax)
	addl	$16, %eax
	cmpl	%ebx, %edx
	jl	L41
L126:
	incl	%ecx
	cmpl	%ebx, %ecx
	jle	L152
L124:
	movl	-348(%ebp), %eax
	fld1
	xorl	%edx, %edx
	movl	%edx, -208(%ebp)
	xorl	%edi, %edi
	leal	-88(%ebp), %ecx
	fstpl	-216(%ebp)
	sall	$4, %eax
	addl	%eax, %ecx
	movl	%edi, -204(%ebp)
	fldl	8(%ecx)
	fldl	(%ecx)
	fld	%st(1)
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(2), %st
	faddp	%st, %st(1)
	fld	%st(2)
	fmul	%st(4), %st
	fadd	%st(2), %st
	fxch	%st(2)
	fmulp	%st, %st(4)
	fdivr	%st, %st(1)
	fxch	%st(3)
	fsubp	%st, %st(2)
	fstl	-216(%ebp)
	fstl	-200(%ebp)
	fxch	%st(1)
	fdivp	%st, %st(2)
	fstl	-120(%ebp)
	fxch	%st(1)
	fstl	-208(%ebp)
	fstl	-192(%ebp)
	fstl	-112(%ebp)
	fxch	%st(1)
	fstl	-104(%ebp)
	fxch	%st(1)
	fstpl	-96(%ebp)
	fldl	(%ecx)
	fldl	_rminshift
	fxch	%st(1)
	fabs
	fxch	%st(1)
	fucom	%st(1)
	fnstsw	%ax
	fstp	%st(1)
	sahf
	jbe	L147
	fldl	8(%ecx)
	fabs
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	ja	L143
L150:
	fstp	%st(0)
	fldl	(%esi)
	movl	$1, %ecx
	movl	72(%ebp), %edi
	fldl	-104(%ebp)
	fld	%st(1)
	fstl	-312(%ebp)
	fxch	%st(2)
	fmul	%st(1), %st
	fldl	8(%esi)
	fldl	-96(%ebp)
	fld	%st(1)
	fmul	%st(1), %st
	fxch	%st(5)
	fmulp	%st, %st(1)
	fxch	%st(3)
	fmulp	%st, %st(1)
	fxch	%st(1)
	fsubp	%st, %st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-312(%ebp)
	fstl	-296(%ebp)
	fxch	%st(1)
	fstl	-288(%ebp)
	fstpl	-304(%ebp)
	faddl	(%edi)
	fstpl	(%edi)
	fldl	-288(%ebp)
	faddl	8(%edi)
	fstpl	8(%edi)
	jmp	L140
	.p2align 4,,7
L144:
	fldl	-96(%ebp)
	movl	%ecx, %eax
	incl	%ecx
	fldl	-104(%ebp)
	movl	72(%ebp), %edi
	sall	$4, %eax
	fldl	-112(%ebp)
	fld	%st(1)
	fld	%st(3)
	fldl	-120(%ebp)
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(2)
	leal	(%eax,%edi), %edx
	fmul	%st(1), %st
	fxch	%st(1)
	addl	%esi, %eax
	fmulp	%st, %st(5)
	fxch	%st(3)
	fmulp	%st, %st(2)
	fsubrp	%st, %st(2)
	faddp	%st, %st(2)
	fstl	-104(%ebp)
	fxch	%st(1)
	fstl	-96(%ebp)
	fldl	(%eax)
	fld	%st(0)
	fmul	%st(3), %st
	fxch	%st(1)
	fstl	-344(%ebp)
	fmul	%st(2), %st
	fldl	8(%eax)
	fld	%st(0)
	fmulp	%st, %st(4)
	fmulp	%st, %st(4)
	fxch	%st(1)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstl	-344(%ebp)
	fstl	-328(%ebp)
	fxch	%st(1)
	fstl	-336(%ebp)
	fstpl	-320(%ebp)
	faddl	(%edx)
	fstpl	(%edx)
	fldl	-320(%ebp)
	faddl	8(%edx)
	fstpl	8(%edx)
L140:
	cmpl	%ebx, %ecx
	jle	L144
L96:
	incl	-348(%ebp)
	movl	%ebx, %ecx
	sall	$4, %ecx
	movl	72(%ebp), %edx
	cmpl	$3, -348(%ebp)
	leal	16(%ecx,%edx), %esi
	movl	%esi, 72(%ebp)
	jg	L145
	fldz
	jmp	L118
L147:
	fstp	%st(0)
	jmp	L150
L145:
	addl	$380, %esp
	popl	%ebx
	popl	%esi
	popl	%edi
	popl	%ebp
	ret
L143:
	fstpl	(%esp)
	call	___fpclassify
	testb	$1, %ah
	je	L146
L53:
	movl	76(%ebp), %ecx
	movl	72(%ebp), %edi
	fldl	(%ecx)
	faddl	(%edi)
	fstpl	(%edi)
	fldl	8(%ecx)
	faddl	8(%edi)
	fstpl	8(%edi)
L67:
	movl	$1, %edx
	movl	_pshift, %ebx
	movl	%edx, -352(%ebp)
	cmpl	%ebx, -352(%ebp)
	jg	L96
	movl	72(%ebp), %esi
	movl	76(%ebp), %edi
	addl	$16, %esi
	addl	$16, %edi
	jmp	L95
L77:
	fldl	(%edi)
	faddl	(%esi)
	fstpl	(%esi)
	fldl	8(%edi)
L139:
	faddl	8(%esi)
	addl	$16, %edi
	movl	_pshift, %ebx
	incl	-352(%ebp)
	fstpl	8(%esi)
	addl	$16, %esi
	cmpl	%ebx, -352(%ebp)
	jg	L96
L95:
	fldl	-104(%ebp)
	fldl	-96(%ebp)
	fldl	-120(%ebp)
	fld	%st(2)
	fld	%st(2)
	fldl	-112(%ebp)
	fxch	%st(2)
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(2), %st
	fxch	%st(5)
	fmulp	%st, %st(2)
	fxch	%st(2)
	fmulp	%st, %st(3)
	fxch	%st(1)
	fsubp	%st, %st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-104(%ebp)
	fstpl	(%esp)
	fstpl	-96(%ebp)
	call	___fpclassify
	testb	$1, %ah
	jne	L77
	fldl	-96(%ebp)
	fstpl	(%esp)
	call	___fpclassify
	testb	$1, %ah
	jne	L77
	fldl	-104(%ebp)
	movl	-352(%ebp), %eax
	movl	_wksp, %ebx
	fldl	-96(%ebp)
	sall	$4, %eax
	addl	%ebx, %eax
	fldl	(%eax)
	fld	%st(0)
	fmul	%st(3), %st
	fxch	%st(1)
	fstl	-280(%ebp)
	fmul	%st(2), %st
	fldl	8(%eax)
	fld	%st(0)
	fmulp	%st, %st(4)
	fmulp	%st, %st(4)
	fxch	%st(1)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstl	-280(%ebp)
	fstl	-264(%ebp)
	fxch	%st(1)
	fstl	-256(%ebp)
	fstpl	-272(%ebp)
	faddl	(%esi)
	fstpl	(%esi)
	fldl	-256(%ebp)
	jmp	L139
L146:
	fldl	-96(%ebp)
	fstpl	(%esp)
	call	___fpclassify
	testb	$1, %ah
	jne	L53
	fldl	-104(%ebp)
	movl	_wksp, %esi
	movl	72(%ebp), %ebx
	fldl	-96(%ebp)
	fldl	(%esi)
	fld	%st(0)
	fmul	%st(3), %st
	fxch	%st(1)
	fstl	-248(%ebp)
	fmul	%st(2), %st
	fldl	8(%esi)
	fld	%st(0)
	fmulp	%st, %st(4)
	fmulp	%st, %st(4)
	fxch	%st(1)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstl	-248(%ebp)
	fstl	-232(%ebp)
	fxch	%st(1)
	fstl	-224(%ebp)
	fstpl	-240(%ebp)
	faddl	(%ebx)
	fstpl	(%ebx)
	fldl	-224(%ebp)
	faddl	8(%ebx)
	fstpl	8(%ebx)
	jmp	L67
	.align 2
	.p2align 4,,15
.globl __Z9shift_m2mSt7complexIdES0_S0_S0_PS0_PKS0_
	.def	__Z9shift_m2mSt7complexIdES0_S0_S0_PS0_PKS0_;	.scl	2;	.type	32;	.endef
__Z9shift_m2mSt7complexIdES0_S0_S0_PS0_PKS0_:
	pushl	%ebp
	movl	$3, %edx
	movl	%esp, %ebp
	pushl	%edi
	leal	-88(%ebp), %eax
	fldz
	pushl	%esi
	pushl	%ebx
	subl	$524, %esp
	movl	76(%ebp), %edi
	cmpl	$-1, %edx
	je	L155
L326:
	movl	$0, (%eax)
	decl	%edx
	movl	$0, 4(%eax)
	movl	$0, 8(%eax)
	movl	$0, 12(%eax)
	addl	$16, %eax
	cmpl	$-1, %edx
	jne	L326
L155:
	fldl	8(%ebp)
	xorl	%eax, %eax
	movl	$0, -104(%ebp)
	movl	_wksp, %esi
	movl	$0, -100(%ebp)
	movl	$0, -96(%ebp)
	fstpl	-88(%ebp)
	fldl	16(%ebp)
	movl	$0, -92(%ebp)
	movl	$0, -120(%ebp)
	movl	$0, -116(%ebp)
	fstpl	-80(%ebp)
	fldl	24(%ebp)
	movl	$0, -112(%ebp)
	movl	$0, -108(%ebp)
	movl	%eax, -476(%ebp)
	fstpl	-72(%ebp)
	fldl	32(%ebp)
	fstpl	-64(%ebp)
	fldl	40(%ebp)
	fstpl	-56(%ebp)
	fldl	48(%ebp)
	fstpl	-48(%ebp)
	fldl	56(%ebp)
	fstpl	-40(%ebp)
	fldl	64(%ebp)
	fstpl	-32(%ebp)
L291:
	fldl	_rminshift
	movl	-476(%ebp), %edx
	leal	-88(%ebp), %ecx
	sall	$4, %edx
	addl	%edx, %ecx
	fldl	(%ecx)
	fabs
	fxch	%st(1)
	fucom	%st(1)
	fnstsw	%ax
	fstp	%st(1)
	sahf
	jbe	L321
	fldl	8(%ecx)
	fabs
	fxch	%st(1)
	fucompp
	fnstsw	%ax
	sahf
	ja	L316
L164:
	movl	-476(%ebp), %eax
	fld1
	xorl	%ecx, %ecx
	movl	%ecx, -240(%ebp)
	xorl	%ebx, %ebx
	leal	-88(%ebp), %edx
	movl	%ebx, -236(%ebp)
	movl	_pshift, %ebx
	sall	$4, %eax
	fstpl	-248(%ebp)
	movl	$2, %ecx
	leal	(%edx,%eax), %eax
	fldl	8(%eax)
	fldl	(%eax)
	fld	%st(1)
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(2), %st
	faddp	%st, %st(1)
	fld	%st(2)
	fmul	%st(4), %st
	fadd	%st(2), %st
	fxch	%st(2)
	fmulp	%st, %st(4)
	fdivr	%st, %st(1)
	fxch	%st(3)
	fsubp	%st, %st(2)
	fstl	-104(%ebp)
	fstl	-248(%ebp)
	fxch	%st(1)
	fdivp	%st, %st(2)
	fstl	-232(%ebp)
	fxch	%st(1)
	fstl	-96(%ebp)
	fstl	-240(%ebp)
	fstl	-224(%ebp)
	fxch	%st(1)
	fstpl	-120(%ebp)
	fstpl	-112(%ebp)
	fldl	(%edi)
	fstpl	(%esi)
	fldl	8(%edi)
	fstpl	8(%esi)
	fldl	16(%edi)
	fldl	-104(%ebp)
	fldl	-96(%ebp)
	fld	%st(2)
	fmul	%st(2), %st
	fxch	%st(3)
	fstl	-280(%ebp)
	fmul	%st(1), %st
	fldl	24(%edi)
	fld	%st(0)
	fmulp	%st, %st(3)
	fmulp	%st, %st(3)
	fxch	%st(3)
	fsubp	%st, %st(1)
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-280(%ebp)
	fstl	-264(%ebp)
	fxch	%st(1)
	fstl	-256(%ebp)
	fstpl	-272(%ebp)
	fstpl	16(%esi)
	fldl	-256(%ebp)
	fstpl	24(%esi)
	jmp	L313
	.p2align 4,,7
L317:
	fldl	-96(%ebp)
	movl	%ecx, %edx
	sall	$4, %edx
	fldl	-104(%ebp)
	incl	%ecx
	fld	%st(1)
	fldl	-112(%ebp)
	fld	%st(2)
	leal	(%edx,%edi), %eax
	fldl	-120(%ebp)
	fxch	%st(3)
	fmul	%st(2), %st
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(4)
	fmulp	%st, %st(2)
	fxch	%st(2)
	fmulp	%st, %st(4)
	fxch	%st(2)
	fsubp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(2)
	fstl	-104(%ebp)
	fxch	%st(1)
	fstl	-96(%ebp)
	fldl	(%eax)
	fld	%st(0)
	fmul	%st(3), %st
	fxch	%st(1)
	fstl	-312(%ebp)
	fmul	%st(2), %st
	fldl	8(%eax)
	fld	%st(0)
	fmulp	%st, %st(4)
	fmulp	%st, %st(4)
	fxch	%st(1)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstl	-312(%ebp)
	fstl	-296(%ebp)
	fxch	%st(1)
	fstl	-304(%ebp)
	fstpl	-288(%ebp)
	fstpl	(%edx,%esi)
	fldl	-288(%ebp)
	fstpl	8(%edx,%esi)
L313:
	cmpl	%ebx, %ecx
	jle	L317
L213:
	cmpl	$1, %ebx
	movl	%ebx, %ecx
	jle	L298
	.p2align 4,,15
L327:
	cmpl	%ebx, %ecx
	movl	%ecx, %edx
	jg	L300
	movl	%ecx, %eax
	sall	$4, %eax
	addl	%esi, %eax
	.p2align 4,,15
L243:
	fldl	-16(%eax)
	incl	%edx
	faddl	(%eax)
	fstpl	(%eax)
	fldl	-8(%eax)
	faddl	8(%eax)
	fstpl	8(%eax)
	addl	$16, %eax
	cmpl	%ebx, %edx
	jle	L243
L300:
	decl	%ecx
	cmpl	$1, %ecx
	jg	L327
L298:
	movl	-476(%ebp), %edx
	leal	-24(%ebp), %ecx
	sall	$4, %edx
	leal	(%edx,%ecx), %eax
	movl	_potshift, %ecx
	fldl	-64(%eax)
	testl	%ecx, %ecx
	fstl	-104(%ebp)
	fldl	-56(%eax)
	movl	72(%ebp), %eax
	fstl	-96(%ebp)
	fldl	(%edi)
	faddl	(%eax)
	fstpl	(%eax)
	fldl	8(%edi)
	faddl	8(%eax)
	fstpl	8(%eax)
	jne	L248
	fstp	%st(0)
	fstp	%st(0)
	movl	$1, %ecx
	cmpl	%ebx, %ecx
	jg	L269
	leal	-88(%ebp), %eax
	addl	%edx, %eax
	fldz
	movl	%eax, -484(%ebp)
	fld1
L268:
	pushl	%ecx
	movl	72(%ebp), %eax
	movl	%ecx, %edx
	fildl	(%esp)
	incl	%ecx
	addl	$4, %esp
	fldl	(%edi)
	sall	$4, %edx
	movl	%edx, -508(%ebp)
	fld	%st(0)
	addl	%eax, %edx
	fstl	-376(%ebp)
	movl	-508(%ebp), %eax
	fmul	%st(4), %st
	fldl	8(%edi)
	addl	%esi, %eax
	cmpl	%ebx, %ecx
	fldl	-96(%ebp)
	fxch	%st(4)
	fdivr	%st(5), %st
	fld	%st(1)
	fxch	%st(4)
	fmul	%st(1), %st
	fxch	%st(2)
	fmulp	%st, %st(1)
	fxch	%st(3)
	fmul	%st(6), %st
	fxch	%st(2)
	faddp	%st, %st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
	fstl	-368(%ebp)
	fstl	-352(%ebp)
	fxch	%st(1)
	fstl	-376(%ebp)
	fstl	-360(%ebp)
	fldl	(%eax)
	fstl	-392(%ebp)
	fsubp	%st, %st(1)
	fldl	8(%eax)
	fldl	-104(%ebp)
	fxch	%st(1)
	fsubp	%st, %st(3)
	fld	%st(1)
	fstl	-392(%ebp)
	fld	%st(3)
	fxch	%st(3)
	fmul	%st(2), %st
	fxch	%st(4)
	fstl	-384(%ebp)
	fxch	%st(3)
	fmul	%st(5), %st
	fxch	%st(1)
	fstl	-344(%ebp)
	fmul	%st(5), %st
	fxch	%st(4)
	fsubp	%st, %st(1)
	fxch	%st(2)
	fstl	-336(%ebp)
	fmul	%st(1), %st
	fld	%st(1)
	fxch	%st(3)
	fstl	-408(%ebp)
	fstl	-328(%ebp)
	fxch	%st(4)
	faddp	%st, %st(1)
	fstl	-400(%ebp)
	fstl	-320(%ebp)
	fxch	%st(3)
	faddl	(%edx)
	fxch	%st(3)
	faddl	8(%edx)
	fxch	%st(3)
	fstpl	(%edx)
	fxch	%st(2)
	fstpl	8(%edx)
	movl	-484(%ebp), %edx
	fld	%st(2)
	fldl	(%edx)
	fldl	8(%edx)
	fxch	%st(3)
	fmul	%st(1), %st
	fxch	%st(2)
	fmul	%st(3), %st
	fxch	%st(4)
	fmulp	%st, %st(3)
	fmulp	%st, %st(4)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstpl	-104(%ebp)
	fstpl	-96(%ebp)
	jle	L268
L322:
	fstp	%st(0)
	fstp	%st(0)
L269:
	incl	-476(%ebp)
	sall	$4, %ebx
	leal	16(%ebx,%edi), %edi
	cmpl	$3, -476(%ebp)
	jg	L319
	fldz
	jmp	L291
L248:
	fldl	16(%esi)
	movl	72(%ebp), %ecx
	movl	72(%ebp), %eax
	fld	%st(0)
	fmul	%st(3), %st
	fxch	%st(1)
	addl	$16, %eax
	fstl	-440(%ebp)
	fmul	%st(2), %st
	fldl	24(%esi)
	fld	%st(0)
	fmul	%st(4), %st
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(2)
	fstl	-440(%ebp)
	fstl	-424(%ebp)
	fxch	%st(1)
	fstl	-432(%ebp)
	fstl	-416(%ebp)
	fxch	%st(1)
	faddl	16(%ecx)
	fstpl	16(%ecx)
	movl	$2, %ecx
	cmpl	%ebx, %ecx
	faddl	8(%eax)
	fstpl	8(%eax)
	jg	L322
	leal	-88(%ebp), %eax
	addl	%edx, %eax
	movl	%eax, -488(%ebp)
	jmp	L290
L323:
	fxch	%st(1)
L290:
	movl	-488(%ebp), %edx
	fld	%st(1)
	movl	%ecx, %eax
	fld	%st(1)
	incl	%ecx
	sall	$4, %eax
	fldl	(%edx)
	fldl	8(%edx)
	fxch	%st(3)
	movl	%eax, %edx
	fmul	%st(1), %st
	fxch	%st(2)
	movl	%eax, -508(%ebp)
	addl	%esi, %eax
	addl	72(%ebp), %edx
	fmul	%st(3), %st
	fxch	%st(3)
	cmpl	%ebx, %ecx
	fmulp	%st, %st(5)
	fmulp	%st, %st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(2)
	fld	%st(0)
	fstpl	-104(%ebp)
	fld	%st(1)
	fstpl	-96(%ebp)
	fldl	(%eax)
	fld	%st(0)
	fmul	%st(2), %st
	fxch	%st(1)
	fstl	-472(%ebp)
	fmul	%st(3), %st
	fldl	8(%eax)
	fld	%st(0)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(2)
	fstl	-472(%ebp)
	fstl	-456(%ebp)
	fxch	%st(1)
	fstl	-464(%ebp)
	fstl	-448(%ebp)
	fxch	%st(1)
	faddl	(%edx)
	fxch	%st(1)
	faddl	8(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fstpl	8(%edx)
	jle	L323
	fstp	%st(0)
	fstp	%st(0)
	jmp	L269
L321:
	fstp	%st(0)
	jmp	L164
L319:
	addl	$524, %esp
	popl	%ebx
	popl	%esi
	popl	%edi
	popl	%ebp
	ret
L316:
	fld1
	xorl	%edx, %edx
	xorl	%ebx, %ebx
	fstpl	-152(%ebp)
	movl	%edx, -144(%ebp)
	movl	%ebx, -140(%ebp)
	fldl	(%ecx)
	fldl	8(%ecx)
	fld	%st(1)
	fld	%st(1)
	fmul	%st(2), %st
	fxch	%st(1)
	fmul	%st(3), %st
	faddp	%st, %st(1)
	fld	%st(1)
	fmul	%st(4), %st
	fadd	%st(3), %st
	fxch	%st(3)
	fmulp	%st, %st(4)
	fdivr	%st, %st(2)
	fxch	%st(3)
	fsubp	%st, %st(1)
	fxch	%st(1)
	fstl	-104(%ebp)
	fstl	-152(%ebp)
	fxch	%st(1)
	fdivp	%st, %st(2)
	fstl	-136(%ebp)
	fxch	%st(1)
	fstl	-144(%ebp)
	fstl	-128(%ebp)
	fxch	%st(1)
	fstpl	-120(%ebp)
	fstl	-112(%ebp)
	fstpl	-96(%ebp)
	fldl	(%edi)
	fstpl	(%esi)
	fldl	8(%edi)
	fstpl	8(%esi)
	fldl	-104(%ebp)
	fstpl	(%esp)
	call	___fpclassify
	testb	$1, %ah
	je	L320
L174:
	movl	_wksp, %esi
	fldz
	movl	72(%ebp), %eax
	movl	72(%ebp), %ecx
	fstl	16(%esi)
	fstpl	24(%esi)
	addl	$16, %ecx
	fldl	16(%edi)
	faddl	16(%eax)
	fstpl	16(%eax)
	fldl	24(%edi)
	faddl	8(%ecx)
	fstpl	8(%ecx)
L185:
	movl	$2, %ebx
	movl	%ebx, -480(%ebp)
	movl	_pshift, %ebx
	cmpl	%ebx, -480(%ebp)
	jg	L213
	leal	32(%edi), %esi
	movl	%esi, -492(%ebp)
	jmp	L212
L196:
	movl	-480(%ebp), %eax
	movl	_wksp, %esi
	movl	-492(%ebp), %ebx
	sall	$4, %eax
	movl	72(%ebp), %ecx
	leal	(%eax,%esi), %edx
	movl	$0, (%edx)
	addl	%ecx, %eax
	movl	$0, 4(%edx)
	movl	$0, 8(%edx)
	movl	$0, 12(%edx)
	fldl	(%ebx)
	faddl	(%eax)
	fstpl	(%eax)
	fldl	8(%ebx)
	faddl	8(%eax)
	fstpl	8(%eax)
L192:
	incl	-480(%ebp)
	movl	_pshift, %ebx
	addl	$16, -492(%ebp)
	cmpl	%ebx, -480(%ebp)
	jg	L213
L212:
	fldl	-104(%ebp)
	fldl	-96(%ebp)
	fldl	-120(%ebp)
	fld	%st(2)
	fld	%st(2)
	fldl	-112(%ebp)
	fxch	%st(2)
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(2), %st
	fxch	%st(5)
	fmulp	%st, %st(2)
	fxch	%st(2)
	fmulp	%st, %st(3)
	fxch	%st(1)
	fsubp	%st, %st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-104(%ebp)
	fstpl	(%esp)
	fstpl	-96(%ebp)
	call	___fpclassify
	testb	$1, %ah
	jne	L196
	fldl	-96(%ebp)
	fstpl	(%esp)
	call	___fpclassify
	testb	$1, %ah
	jne	L196
	fldl	-104(%ebp)
	movl	-492(%ebp), %edx
	movl	-480(%ebp), %ecx
	fldl	-96(%ebp)
	movl	_wksp, %esi
	sall	$4, %ecx
	fldl	(%edx)
	fld	%st(0)
	fmul	%st(3), %st
	fxch	%st(1)
	fstl	-216(%ebp)
	fmul	%st(2), %st
	fldl	8(%edx)
	fld	%st(0)
	fmulp	%st, %st(4)
	fmulp	%st, %st(4)
	fxch	%st(1)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstl	-216(%ebp)
	fstl	-200(%ebp)
	fxch	%st(1)
	fstl	-192(%ebp)
	fstpl	-208(%ebp)
	fstpl	(%ecx,%esi)
	fldl	-192(%ebp)
	fstpl	8(%ecx,%esi)
	jmp	L192
L320:
	fldl	-96(%ebp)
	fstpl	(%esp)
	call	___fpclassify
	testb	$1, %ah
	jne	L174
	fldl	16(%edi)
	movl	_wksp, %esi
	fldl	-104(%ebp)
	fld	%st(1)
	fstl	-184(%ebp)
	fxch	%st(2)
	fmul	%st(1), %st
	fldl	24(%edi)
	fldl	-96(%ebp)
	fld	%st(1)
	fmul	%st(1), %st
	fxch	%st(5)
	fmulp	%st, %st(1)
	fxch	%st(3)
	fmulp	%st, %st(1)
	fxch	%st(1)
	fsubp	%st, %st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-184(%ebp)
	fstl	-168(%ebp)
	fxch	%st(1)
	fstl	-160(%ebp)
	fstpl	-176(%ebp)
	fstpl	16(%esi)
	fldl	-160(%ebp)
	fstpl	24(%esi)
	jmp	L185
	.section .rdata,"dr"
	.align 8
LC18:
	.long	1413754136
	.long	1074340347
	.text
	.align 2
	.p2align 4,,15
.globl __Z10shift_m2psSt7complexIdEPKvjPKiiiPKS0_PS0_S6_S7_
	.def	__Z10shift_m2psSt7complexIdEPKvjPKiiiPKS0_PS0_S6_S7_;	.scl	2;	.type	32;	.endef
__Z10shift_m2psSt7complexIdEPKvjPKiiiPKS0_PS0_S6_S7_:
	pushl	%ebp
	movl	%esp, %ebp
	pushl	%edi
	pushl	%esi
	pushl	%ebx
	subl	$1420, %esp
	movl	36(%ebp), %eax
	fldl	8(%ebp)
	movl	40(%ebp), %edx
	fldl	16(%ebp)
	movl	%eax, -1284(%ebp)
	cmpl	%edx, %eax
	jge	L649
	movl	_pshift, %ecx
	movl	%ecx, -1404(%ebp)
L581:
	movl	-1284(%ebp), %eax
	fldz
	movl	32(%ebp), %ebx
	movl	28(%ebp), %edi
	fld1
	movl	24(%ebp), %esi
	movl	(%ebx,%eax,4), %edx
	fldz
	movl	_rpow, %ebx
	movl	%edx, -1400(%ebp)
	movl	%edx, %ecx
	imull	%edi, %ecx
	movl	%ebx, -1300(%ebp)
	movl	-1404(%ebp), %edi
	movl	_wksp, %ebx
	movl	-1400(%ebp), %eax
	addl	%esi, %ecx
	sall	$4, %edi
	fldl	8(%ecx)
	fxch	%st(5)
	leal	16(%edi,%ebx), %esi
	addl	$16, %edi
	fsubrl	(%ecx)
	fxch	%st(5)
	movl	$1, %ecx
	fsubp	%st, %st(4)
	imull	%eax, %edi
	movl	52(%ebp), %edx
	fld	%st(4)
	fstl	-64(%ebp)
	fld	%st(4)
	addl	%edx, %edi
	fstl	-56(%ebp)
	fmul	%st, %st(5)
	fxch	%st(6)
	movl	%edi, -1412(%ebp)
	fmul	%st(1), %st
	fxch	%st(4)
	fmul	%st(1), %st
	fxch	%st(1)
	fstl	-48(%ebp)
	fxch	%st(4)
	faddp	%st, %st(5)
	fldz
	fxch	%st(6)
	fstl	-40(%ebp)
	fmul	%st, %st(6)
	fsubrp	%st, %st(1)
	fxch	%st(5)
	faddp	%st, %st(3)
	fxch	%st(4)
	fdiv	%st(3), %st
	fxch	%st(2)
	fdivp	%st, %st(3)
	fxch	%st(1)
	fstl	-72(%ebp)
	fxch	%st(2)
	fstl	-1296(%ebp)
	fstl	-80(%ebp)
	fstpl	-32(%ebp)
	fxch	%st(1)
	fstl	-24(%ebp)
	jmp	L630
	.p2align 4,,7
L639:
	fldl	-1296(%ebp)
	movl	%ecx, %edx
	fld	%st(3)
	fmull	-1296(%ebp)
	fxch	%st(1)
	movl	-1300(%ebp), %edi
	sall	$4, %edx
	fmul	%st(3), %st
	fxch	%st(4)
	fmul	%st(2), %st
	fxch	%st(3)
	leal	(%edx,%edi), %eax
	fmul	%st(2), %st
	fxch	%st(4)
	fsubp	%st, %st(3)
	faddp	%st, %st(3)
	fxch	%st(1)
	fstpl	(%edx,%edi)
	fxch	%st(1)
	movl	-1412(%ebp), %edi
	fstl	8(%eax)
	fldl	(%eax)
	leal	(%edx,%edi), %eax
	movl	44(%ebp), %edi
	fldl	8(%eax)
	fldl	(%eax)
	leal	(%edx,%ebx), %eax
	fld	%st(1)
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(3)
	fmul	%st(4), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(2)
	fmul	%st(5), %st
	fxch	%st(2)
	fsubp	%st, %st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
	fchs
	fxch	%st(1)
	fchs
	fxch	%st(1)
	fstpl	-16(%eax)
	fstpl	-8(%eax)
	leal	(%edx,%edi), %eax
	addl	%esi, %edx
	fldl	8(%eax)
	movl	-1412(%ebp), %edi
	fldl	(%eax)
	fld	%st(1)
	movl	-1300(%ebp), %eax
	fmul	%st(3), %st
	fld	%st(1)
	fmul	%st(5), %st
	fxch	%st(3)
	fmul	%st(5), %st
	fxch	%st(2)
	fmul	%st(4), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	-1296(%ebp)
	fxch	%st(3)
	fsubp	%st, %st(2)
	fld	%st(4)
	fmull	-1296(%ebp)
	fxch	%st(3)
	fmul	%st(4), %st
	fxch	%st(5)
	fmul	%st(6), %st
	fxch	%st(2)
	fstpl	-16(%edx)
	fxch	%st(3)
	fmul	%st(5), %st
	fxch	%st(4)
	fsubp	%st, %st(1)
	fxch	%st(2)
	fstpl	-8(%edx)
	leal	1(%ecx), %edx
	faddp	%st, %st(2)
	sall	$4, %edx
	addl	$2, %ecx
	fstpl	(%edx,%eax)
	leal	(%edx,%eax), %eax
	fld	%st(0)
	fldl	(%eax)
	fxch	%st(2)
	fstl	8(%eax)
	leal	(%edx,%edi), %eax
	movl	44(%ebp), %edi
	fldl	(%eax)
	fldl	8(%eax)
	leal	(%edx,%ebx), %eax
	fld	%st(1)
	fld	%st(1)
	fmul	%st(4), %st
	fxch	%st(2)
	fmul	%st(6), %st
	fxch	%st(1)
	fmul	%st(6), %st
	fxch	%st(3)
	fmul	%st(4), %st
	fxch	%st(3)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstpl	-16(%eax)
	fstpl	-8(%eax)
	leal	(%edx,%edi), %eax
	addl	%esi, %edx
	fldl	(%eax)
	fldl	8(%eax)
	fld	%st(1)
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(6), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(3)
	fmulp	%st, %st(4)
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(1)
	fsubp	%st, %st(2)
	faddp	%st, %st(2)
	fstpl	-16(%edx)
	fstpl	-8(%edx)
	fxch	%st(2)
L630:
	cmpl	-1404(%ebp), %ecx
	jl	L639
	testb	$1, -1404(%ebp)
	je	L650
	fldl	-1296(%ebp)
	fld	%st(3)
	movl	-1404(%ebp), %edx
	fmull	-1296(%ebp)
	movl	-1300(%ebp), %eax
	fld	%st(3)
	fxch	%st(2)
	fmulp	%st, %st(4)
	fxch	%st(4)
	movl	-1412(%ebp), %ecx
	sall	$4, %edx
	fmul	%st(2), %st
	fxch	%st(1)
	leal	(%edx,%eax), %edi
	fmulp	%st, %st(2)
	fxch	%st(3)
	fstl	-1296(%ebp)
	fxch	%st(2)
	fsubp	%st, %st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstpl	(%edx,%eax)
	leal	(%edx,%ecx), %eax
	movl	44(%ebp), %ecx
	fstl	8(%edi)
	fldl	(%eax)
	fldl	8(%eax)
	leal	(%edx,%ecx), %eax
	fld	%st(1)
	fldl	(%edi)
	fld	%st(2)
	fmul	%st(5), %st
	fxch	%st(3)
	leal	(%edx,%ebx), %edi
	addl	%esi, %edx
	fmul	%st(1), %st
	fxch	%st(2)
	fmul	%st(1), %st
	fxch	%st(4)
	fmul	%st(5), %st
	fxch	%st(4)
	fsubp	%st, %st(3)
	fxch	%st(3)
	faddp	%st, %st(1)
	fxch	%st(1)
	fchs
	fxch	%st(1)
	fchs
	fxch	%st(1)
	fstpl	-16(%edi)
	fstpl	-8(%edi)
	fldl	(%eax)
	fldl	8(%eax)
	fld	%st(1)
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(3)
	fmulp	%st, %st(5)
	fxch	%st(1)
	fmulp	%st, %st(3)
	fsubp	%st, %st(1)
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstpl	-16(%edx)
	fstpl	-8(%edx)
L362:
	movl	-1404(%ebp), %ecx
	movl	$2, %edi
	movl	-1404(%ebp), %edx
	movl	%edi, -1304(%ebp)
	sall	$4, %ecx
	subl	$3, %edx
	leal	(%ecx,%ebx), %eax
	addl	%esi, %ecx
	cmpl	$2, %edx
	movl	$0, 8(%eax)
	movl	$0, 12(%eax)
	movl	$0, (%eax)
	movl	$0, 4(%eax)
	movl	$0, 8(%ecx)
	movl	$0, 12(%ecx)
	movl	$0, (%ecx)
	movl	$0, 4(%ecx)
	jl	L587
	movl	-1404(%ebp), %edi
	subl	$7, %edi
	.p2align 4,,15
L381:
	movl	-1404(%ebp), %ecx
	movl	-1304(%ebp), %edx
	subl	%edx, %ecx
	jmp	L631
	.p2align 4,,7
L640:
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$8, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	32(%eax)
	fldl	(%eax)
	fldl	-16(%eax)
	fldl	-32(%eax)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fldl	-48(%eax)
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(1), %st
	fldl	64(%eax)
	fxch	%st(4)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(3)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(3), %st
	fxch	%st(6)
	fstpl	-48(%eax)
	fxch	%st(5)
	fstpl	-32(%eax)
	fldl	48(%eax)
	fldl	96(%eax)
	fxch	%st(5)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(6), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	80(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fldl	112(%eax)
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(4)
	faddl	128(%eax)
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(5)
	fadd	%st(2), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(2)
	fadd	%st(6), %st
	fxch	%st(6)
	fstpl	112(%eax)
	fadd	%st(4), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(1)
	fstpl	16(%eax)
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fldl	24(%eax)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(3)
	fstpl	32(%eax)
	fldl	8(%eax)
	fldl	40(%eax)
	fxch	%st(5)
	fstpl	80(%eax)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fldl	-8(%eax)
	fxch	%st(2)
	fstpl	48(%eax)
	fldl	-24(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(6)
	fstpl	96(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(5), %st
	fldl	-40(%eax)
	fxch	%st(4)
	fstpl	64(%eax)
	fldl	72(%eax)
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(3)
	fstpl	-40(%eax)
	fxch	%st(2)
	fstpl	-24(%eax)
	fldl	56(%eax)
	fadd	%st, %st(4)
	fadd	%st(3), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(5), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(4)
	fstpl	-8(%eax)
	fxch	%st(3)
	fstpl	8(%eax)
	fldl	88(%eax)
	fldl	104(%eax)
	fldl	120(%eax)
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(5)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(4)
	faddl	136(%eax)
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(6)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(5)
	fadd	%st(1), %st
	fxch	%st(4)
	fadd	%st(6), %st
	fxch	%st(1)
	fstpl	120(%eax)
	fxch	%st(5)
	fadd	%st(3), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(3)
	fstpl	104(%eax)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(4), %st
	fxch	%st(2)
	fstpl	24(%eax)
	fadd	%st, %st(3)
	fxch	%st(2)
	fstpl	40(%eax)
	fstpl	56(%eax)
	fstpl	88(%eax)
	fstpl	72(%eax)
	fldl	16(%edx)
	fldl	32(%edx)
	fldl	(%edx)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fldl	-48(%edx)
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(1), %st
	fldl	64(%edx)
	fxch	%st(4)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(3)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(3), %st
	fxch	%st(6)
	fstpl	-48(%edx)
	fxch	%st(5)
	fstpl	-32(%edx)
	fldl	48(%edx)
	fldl	96(%edx)
	fxch	%st(5)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(6), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(1)
	fstpl	-16(%edx)
	fstpl	(%edx)
	fldl	80(%edx)
	fldl	112(%edx)
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(4)
	faddl	128(%edx)
	fxch	%st(5)
	fadd	%st(2), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(4)
	fadd	%st(6), %st
	fxch	%st(5)
	fstpl	112(%edx)
	fxch	%st(5)
	fadd	%st(3), %st
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(3)
	fstpl	96(%edx)
	fadd	%st, %st(4)
	fadd	%st(3), %st
	fxch	%st(2)
	fstpl	16(%edx)
	fadd	%st, %st(2)
	fxch	%st(3)
	fstpl	32(%edx)
	fstpl	48(%edx)
	fxch	%st(1)
	fstpl	80(%edx)
	fstpl	64(%edx)
	fldl	24(%edx)
	fldl	40(%edx)
	fldl	8(%edx)
	fldl	-8(%edx)
	fldl	-24(%edx)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fldl	-40(%edx)
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(1), %st
	fldl	72(%edx)
	fxch	%st(4)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(3)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(3), %st
	fxch	%st(6)
	fstpl	-40(%edx)
	fxch	%st(5)
	fstpl	-24(%edx)
	fldl	56(%edx)
	fldl	104(%edx)
	fxch	%st(5)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(6), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(1)
	fstpl	-8(%edx)
	fstpl	8(%edx)
	fldl	88(%edx)
	fldl	120(%edx)
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(4)
	faddl	136(%edx)
	fxch	%st(5)
	fadd	%st(2), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(4)
	fadd	%st(6), %st
	fxch	%st(5)
	fstpl	120(%edx)
	fxch	%st(5)
	fadd	%st(3), %st
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(3)
	fstpl	104(%edx)
	fadd	%st, %st(4)
	fadd	%st(3), %st
	fxch	%st(2)
	fstpl	24(%edx)
	fadd	%st, %st(2)
	fxch	%st(3)
	fstpl	40(%edx)
	fstpl	56(%edx)
	fxch	%st(1)
	fstpl	88(%edx)
	fstpl	72(%edx)
L631:
	cmpl	%ecx, %edi
	jg	L640
	movl	-1404(%ebp), %edx
	subl	$3, %edx
	cmpl	%ecx, %edx
	movl	%edx, -1308(%ebp)
	jle	L378
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$4, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fldl	-16(%eax)
	fldl	-32(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fldl	-48(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fadd	%st, %st(3)
	fxch	%st(2)
	fadd	%st(3), %st
	fstpl	-48(%eax)
	fldl	32(%eax)
	fadd	%st, %st(4)
	fxch	%st(1)
	fadd	%st(4), %st
	fadd	%st, %st(2)
	fxch	%st(3)
	fadd	%st(2), %st
	fstpl	-32(%eax)
	fldl	48(%eax)
	fadd	%st, %st(1)
	faddl	64(%eax)
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fstpl	48(%eax)
	fadd	%st, %st(2)
	fadd	%st(3), %st
	fldl	-24(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%eax)
	fldl	24(%eax)
	fxch	%st(3)
	fstpl	-16(%eax)
	fstpl	(%eax)
	fldl	8(%eax)
	fldl	-8(%eax)
	fxch	%st(4)
	fstpl	32(%eax)
	fadd	%st(2), %st
	fldl	-40(%eax)
	fxch	%st(4)
	fadd	%st(1), %st
	fadd	%st, %st(2)
	fxch	%st(4)
	fadd	%st(2), %st
	fstpl	-40(%eax)
	fldl	40(%eax)
	fadd	%st, %st(3)
	fxch	%st(1)
	fadd	%st(3), %st
	fadd	%st, %st(4)
	fxch	%st(2)
	fadd	%st(4), %st
	fstpl	-24(%eax)
	fldl	56(%eax)
	fadd	%st, %st(1)
	faddl	72(%eax)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fstpl	56(%eax)
	fadd	%st, %st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	40(%eax)
	fadd	%st, %st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	24(%eax)
	fxch	%st(1)
	fstpl	-8(%eax)
	fstpl	8(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fldl	-48(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fadd	%st, %st(3)
	fxch	%st(2)
	fadd	%st(3), %st
	fstpl	-48(%edx)
	fldl	32(%edx)
	fadd	%st, %st(4)
	fxch	%st(1)
	fadd	%st(4), %st
	fadd	%st, %st(2)
	fxch	%st(3)
	fadd	%st(2), %st
	fstpl	-32(%edx)
	fldl	48(%edx)
	fadd	%st, %st(1)
	faddl	64(%edx)
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fstpl	48(%edx)
	fadd	%st, %st(2)
	fadd	%st(3), %st
	fldl	-24(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%edx)
	fldl	24(%edx)
	fxch	%st(3)
	fstpl	-16(%edx)
	fstpl	(%edx)
	fldl	8(%edx)
	fldl	-8(%edx)
	fxch	%st(4)
	fstpl	32(%edx)
	fadd	%st(2), %st
	fldl	-40(%edx)
	fxch	%st(4)
	fadd	%st(1), %st
	fadd	%st, %st(2)
	fxch	%st(4)
	fadd	%st(2), %st
	fstpl	-40(%edx)
	fldl	40(%edx)
	fadd	%st, %st(3)
	fxch	%st(1)
	fadd	%st(3), %st
	fadd	%st, %st(4)
	fxch	%st(2)
	fadd	%st(4), %st
	fstpl	-24(%edx)
	fldl	56(%edx)
	fadd	%st, %st(1)
	faddl	72(%edx)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fstpl	56(%edx)
	fadd	%st, %st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	40(%edx)
	fadd	%st, %st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	24(%edx)
	fxch	%st(1)
	fstpl	-8(%edx)
	fstpl	8(%edx)
L378:
	movl	-1404(%ebp), %edx
	decl	%edx
	cmpl	%ecx, %edx
	jle	L379
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$2, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fldl	-16(%eax)
	fldl	-32(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	faddl	32(%eax)
	fldl	-48(%eax)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	(%eax)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fldl	24(%eax)
	fldl	8(%eax)
	fxch	%st(3)
	fstpl	-16(%eax)
	fldl	-8(%eax)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	40(%eax)
	fxch	%st(2)
	fstpl	-32(%eax)
	fadd	%st, %st(2)
	fldl	-24(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(4)
	fstpl	-48(%eax)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(3), %st
	fldl	-40(%eax)
	fxch	%st(2)
	fstpl	24(%eax)
	fxch	%st(3)
	fstpl	8(%eax)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-8(%eax)
	fstpl	-40(%eax)
	fstpl	-24(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	faddl	32(%edx)
	fldl	-48(%edx)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fldl	24(%edx)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(3)
	fstpl	(%edx)
	fldl	8(%edx)
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	-8(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(3)
	fstpl	-32(%edx)
	fldl	-24(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(4)
	fstpl	-48(%edx)
	fldl	-40(%edx)
	fxch	%st(1)
	fadd	%st(4), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-40(%edx)
	fxch	%st(1)
	faddl	40(%edx)
	fadd	%st, %st(2)
	fstpl	24(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	8(%edx)
	fadd	%st, %st(1)
	fstpl	-8(%edx)
	fstpl	-24(%edx)
L379:
	cmpl	-1404(%ebp), %ecx
	jge	L380
	movl	%ecx, %edx
	sall	$4, %edx
	incl	%ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	faddl	(%eax)
	fstl	(%eax)
	faddl	-16(%eax)
	fstl	-16(%eax)
	faddl	-32(%eax)
	fstl	-32(%eax)
	faddl	-48(%eax)
	fstpl	-48(%eax)
	fldl	24(%eax)
	faddl	8(%eax)
	fstl	8(%eax)
	faddl	-8(%eax)
	fstl	-8(%eax)
	faddl	-24(%eax)
	fstl	-24(%eax)
	faddl	-40(%eax)
	fstpl	-40(%eax)
	fldl	16(%edx)
	faddl	(%edx)
	fstl	(%edx)
	faddl	-16(%edx)
	fstl	-16(%edx)
	faddl	-32(%edx)
	fstl	-32(%edx)
	faddl	-48(%edx)
	fstpl	-48(%edx)
	fldl	24(%edx)
	faddl	8(%edx)
	fstl	8(%edx)
	faddl	-8(%edx)
	fstl	-8(%edx)
	faddl	-24(%edx)
	fstl	-24(%edx)
	faddl	-40(%edx)
	fstpl	-40(%edx)
L380:
	addl	$4, -1304(%ebp)
	leal	-1(%ecx), %edx
	sall	$4, %edx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	movl	-1304(%ebp), %ecx
	fldl	(%eax)
	fldl	-16(%eax)
	fldl	-32(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	cmpl	%ecx, -1308(%ebp)
	fadd	%st, %st(1)
	fadd	%st(3), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(3)
	fldl	24(%eax)
	fxch	%st(1)
	fstpl	-32(%eax)
	fxch	%st(1)
	fstpl	-16(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fldl	8(%eax)
	fldl	-8(%eax)
	fldl	-24(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fadd	%st, %st(1)
	fadd	%st(3), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(3)
	fstpl	-24(%eax)
	fstpl	-8(%eax)
	fstpl	8(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fadd	%st, %st(1)
	fadd	%st(3), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(3)
	fldl	24(%edx)
	fxch	%st(1)
	fstpl	-32(%edx)
	fxch	%st(1)
	fstpl	-16(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fldl	8(%edx)
	fldl	-8(%edx)
	fldl	-24(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fadd	%st, %st(1)
	fadd	%st(3), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(3)
	fstpl	-24(%edx)
	fstpl	-8(%edx)
	fstpl	8(%edx)
	jge	L381
L587:
	movl	-1404(%ebp), %edi
	decl	%edi
	cmpl	-1304(%ebp), %edi
	jl	L382
	movl	-1404(%ebp), %ecx
	movl	-1304(%ebp), %edx
	movl	-1404(%ebp), %edi
	subl	%edx, %ecx
	subl	$7, %edi
	jmp	L632
L641:
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$8, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fldl	-16(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	32(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	(%eax)
	fldl	48(%eax)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	16(%eax)
	fldl	64(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	32(%eax)
	fldl	80(%eax)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	48(%eax)
	fldl	96(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	64(%eax)
	fldl	112(%eax)
	fadd	%st, %st(1)
	faddl	128(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	112(%eax)
	fldl	24(%eax)
	fxch	%st(1)
	fstpl	80(%eax)
	fxch	%st(1)
	fstpl	96(%eax)
	fldl	8(%eax)
	fldl	-8(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-8(%eax)
	fldl	40(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	8(%eax)
	fldl	56(%eax)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	24(%eax)
	fldl	72(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	40(%eax)
	fldl	88(%eax)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	56(%eax)
	fldl	104(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	72(%eax)
	fldl	120(%eax)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	88(%eax)
	fxch	%st(1)
	faddl	136(%eax)
	fadd	%st, %st(1)
	fstpl	120(%eax)
	fstpl	104(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fldl	-16(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	32(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	(%edx)
	fldl	48(%edx)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	16(%edx)
	fldl	64(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	32(%edx)
	fldl	80(%edx)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	48(%edx)
	fldl	96(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	64(%edx)
	fldl	112(%edx)
	fadd	%st, %st(1)
	faddl	128(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	112(%edx)
	fldl	24(%edx)
	fxch	%st(1)
	fstpl	80(%edx)
	fxch	%st(1)
	fstpl	96(%edx)
	fldl	8(%edx)
	fldl	-8(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-8(%edx)
	fldl	40(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	8(%edx)
	fldl	56(%edx)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	24(%edx)
	fldl	72(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	40(%edx)
	fldl	88(%edx)
	fadd	%st, %st(1)
	fxch	%st(2)
	fadd	%st(1), %st
	fstpl	56(%edx)
	fldl	104(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	72(%edx)
	fldl	120(%edx)
	fadd	%st, %st(1)
	faddl	136(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	120(%edx)
	fstpl	88(%edx)
	fstpl	104(%edx)
L632:
	cmpl	%ecx, %edi
	jg	L641
	movl	-1404(%ebp), %edi
	subl	$3, %edi
	cmpl	%ecx, %edi
	jle	L387
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$4, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fldl	-16(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	32(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	(%eax)
	fldl	48(%eax)
	fadd	%st, %st(1)
	faddl	64(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	48(%eax)
	fstpl	16(%eax)
	fldl	8(%eax)
	fldl	24(%eax)
	fxch	%st(2)
	fstpl	32(%eax)
	fldl	-8(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-8(%eax)
	fldl	40(%eax)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	8(%eax)
	fldl	56(%eax)
	fadd	%st, %st(1)
	faddl	72(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	56(%eax)
	fstpl	24(%eax)
	fstpl	40(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fldl	-16(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	32(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	(%edx)
	fldl	48(%edx)
	fadd	%st, %st(1)
	faddl	64(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	48(%edx)
	fstpl	16(%edx)
	fldl	8(%edx)
	fldl	24(%edx)
	fxch	%st(2)
	fstpl	32(%edx)
	fldl	-8(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-8(%edx)
	fldl	40(%edx)
	fadd	%st, %st(2)
	fxch	%st(1)
	fadd	%st(2), %st
	fstpl	8(%edx)
	fldl	56(%edx)
	fadd	%st, %st(1)
	faddl	72(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	56(%edx)
	fstpl	24(%edx)
	fstpl	40(%edx)
L387:
	movl	-1404(%ebp), %edi
	decl	%edi
	cmpl	%ecx, %edi
	jle	L388
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$2, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fldl	-16(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	32(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%eax)
	fldl	24(%eax)
	fxch	%st(2)
	fstpl	-16(%eax)
	fstpl	(%eax)
	fldl	8(%eax)
	fldl	-8(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	40(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	24(%eax)
	fxch	%st(1)
	fstpl	-8(%eax)
	fstpl	8(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fldl	-16(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	32(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%edx)
	fldl	24(%edx)
	fxch	%st(2)
	fstpl	-16(%edx)
	fstpl	(%edx)
	fldl	8(%edx)
	fldl	-8(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	40(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	24(%edx)
	fxch	%st(1)
	fstpl	-8(%edx)
	fstpl	8(%edx)
L388:
	cmpl	-1404(%ebp), %ecx
	jge	L389
	movl	%ecx, %edi
	sall	$4, %edi
	incl	%ecx
	leal	(%edi,%ebx), %edx
	addl	%esi, %edi
	fldl	16(%edx)
	faddl	(%edx)
	fstl	(%edx)
	faddl	-16(%edx)
	fstpl	-16(%edx)
	fldl	24(%edx)
	faddl	8(%edx)
	fstl	8(%edx)
	faddl	-8(%edx)
	fstpl	-8(%edx)
	fldl	16(%edi)
	faddl	(%edi)
	fstl	(%edi)
	faddl	-16(%edi)
	fstpl	-16(%edi)
	fldl	24(%edi)
	faddl	8(%edi)
	fstl	8(%edi)
	faddl	-8(%edi)
	fstpl	-8(%edi)
L389:
	addl	$2, -1304(%ebp)
	leal	-1(%ecx), %eax
	sall	$4, %eax
	leal	(%eax,%ebx), %ecx
	addl	%esi, %eax
	fldl	16(%ecx)
	faddl	(%ecx)
	fstpl	(%ecx)
	fldl	24(%ecx)
	faddl	8(%ecx)
	fstpl	8(%ecx)
	fldl	16(%eax)
	faddl	(%eax)
	fstpl	(%eax)
	fldl	24(%eax)
	faddl	8(%eax)
	fstpl	8(%eax)
L382:
	movl	-1404(%ebp), %ecx
	cmpl	%ecx, -1304(%ebp)
	jg	L390
	movl	-1304(%ebp), %edi
	subl	%edi, %ecx
	movl	-1404(%ebp), %edi
	subl	$7, %edi
	jmp	L633
L642:
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$8, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fadd	%st(1), %st
	fstpl	(%eax)
	fldl	32(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	16(%eax)
	fldl	48(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	32(%eax)
	fldl	64(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	48(%eax)
	fldl	80(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	64(%eax)
	fldl	96(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	80(%eax)
	fldl	112(%eax)
	fadd	%st, %st(1)
	faddl	128(%eax)
	fxch	%st(1)
	fstpl	96(%eax)
	fldl	8(%eax)
	fxch	%st(1)
	fstpl	112(%eax)
	fldl	24(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	8(%eax)
	fldl	40(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	24(%eax)
	fldl	56(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	40(%eax)
	fldl	72(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	56(%eax)
	fldl	88(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	72(%eax)
	fldl	104(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	88(%eax)
	fldl	120(%eax)
	fadd	%st, %st(1)
	faddl	136(%eax)
	fxch	%st(1)
	fstpl	104(%eax)
	fstpl	120(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fadd	%st(1), %st
	fstpl	(%edx)
	fldl	32(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	16(%edx)
	fldl	48(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	32(%edx)
	fldl	64(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	48(%edx)
	fldl	80(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	64(%edx)
	fldl	96(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	80(%edx)
	fldl	112(%edx)
	fadd	%st, %st(1)
	faddl	128(%edx)
	fxch	%st(1)
	fstpl	96(%edx)
	fldl	8(%edx)
	fxch	%st(1)
	fstpl	112(%edx)
	fldl	24(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	8(%edx)
	fldl	40(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	24(%edx)
	fldl	56(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	40(%edx)
	fldl	72(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	56(%edx)
	fldl	88(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	72(%edx)
	fldl	104(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	88(%edx)
	fldl	120(%edx)
	fadd	%st, %st(1)
	faddl	136(%edx)
	fxch	%st(1)
	fstpl	104(%edx)
	fstpl	120(%edx)
L633:
	cmpl	%ecx, %edi
	jg	L642
	movl	-1404(%ebp), %edi
	subl	$3, %edi
	cmpl	%ecx, %edi
	jle	L395
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$4, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fadd	%st(1), %st
	fstpl	(%eax)
	fldl	32(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	16(%eax)
	fldl	48(%eax)
	fadd	%st, %st(1)
	faddl	64(%eax)
	fxch	%st(1)
	fstpl	32(%eax)
	fldl	8(%eax)
	fxch	%st(1)
	fstpl	48(%eax)
	fldl	24(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	8(%eax)
	fldl	40(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	24(%eax)
	fldl	56(%eax)
	fadd	%st, %st(1)
	faddl	72(%eax)
	fxch	%st(1)
	fstpl	40(%eax)
	fstpl	56(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fadd	%st(1), %st
	fstpl	(%edx)
	fldl	32(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	16(%edx)
	fldl	48(%edx)
	fadd	%st, %st(1)
	faddl	64(%edx)
	fxch	%st(1)
	fstpl	32(%edx)
	fldl	8(%edx)
	fxch	%st(1)
	fstpl	48(%edx)
	fldl	24(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	8(%edx)
	fldl	40(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	24(%edx)
	fldl	56(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	40(%edx)
	faddl	72(%edx)
	fstpl	56(%edx)
L395:
	movl	-1404(%ebp), %edi
	decl	%edi
	cmpl	%ecx, %edi
	jle	L396
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$2, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	16(%eax)
	fldl	(%eax)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	32(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fldl	8(%eax)
	fxch	%st(1)
	fstpl	16(%eax)
	fldl	24(%eax)
	fadd	%st, %st(1)
	faddl	40(%eax)
	fxch	%st(1)
	fstpl	8(%eax)
	fstpl	24(%eax)
	fldl	16(%edx)
	fldl	(%edx)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	32(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fldl	8(%edx)
	fxch	%st(1)
	fstpl	16(%edx)
	fldl	24(%edx)
	fadd	%st, %st(1)
	faddl	40(%edx)
	fxch	%st(1)
	fstpl	8(%edx)
	fstpl	24(%edx)
L396:
	cmpl	-1404(%ebp), %ecx
	jge	L390
	movl	%ecx, %edi
	sall	$4, %edi
	leal	(%edi,%ebx), %ecx
	fldl	16(%ecx)
	addl	%esi, %edi
	faddl	(%ecx)
	fstpl	(%ecx)
	fldl	24(%ecx)
	faddl	8(%ecx)
	fstpl	8(%ecx)
	fldl	16(%edi)
	faddl	(%edi)
	fstpl	(%edi)
	fldl	24(%edi)
	faddl	8(%edi)
	fstpl	8(%edi)
L390:
	movl	-1404(%ebp), %edi
	cmpl	$3, %edi
	jle	L595
	movl	%edi, %edx
	subl	$7, %edx
	movl	%edx, -1392(%ebp)
	.p2align 4,,15
L408:
	cmpl	%edi, -1392(%ebp)
	movl	%edi, %ecx
	jmp	L634
	.p2align 4,,7
L643:
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$8, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	-32(%eax)
	fldl	-48(%eax)
	fldl	(%eax)
	fldl	16(%eax)
	fldl	48(%eax)
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fldl	64(%eax)
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(6), %st
	fxch	%st(4)
	faddl	-64(%eax)
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(6), %st
	fxch	%st(6)
	fstpl	-48(%eax)
	fldl	32(%eax)
	fxch	%st(6)
	fstl	-32(%eax)
	fxch	%st(6)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(6)
	fxch	%st(2)
	fadd	%st(3), %st
	fldl	112(%eax)
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(6)
	fstl	-16(%eax)
	faddp	%st, %st(2)
	fadd	%st, %st(4)
	fldl	80(%eax)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(6), %st
	fxch	%st(2)
	fstl	(%eax)
	faddp	%st, %st(6)
	fadd	%st(4), %st
	fldl	8(%eax)
	fxch	%st(5)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(6)
	fstl	16(%eax)
	faddp	%st, %st(2)
	fldl	96(%eax)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(2)
	fstl	32(%eax)
	fxch	%st(3)
	fadd	%st(6), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(6)
	faddp	%st, %st(3)
	fstpl	112(%eax)
	fldl	-40(%eax)
	fxch	%st(3)
	fadd	%st(2), %st
	fldl	-24(%eax)
	fxch	%st(6)
	fstpl	80(%eax)
	fldl	-8(%eax)
	fxch	%st(2)
	fstpl	96(%eax)
	fldl	24(%eax)
	fxch	%st(5)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(6), %st
	fxch	%st(3)
	fstpl	48(%eax)
	fxch	%st(5)
	fadd	%st(3), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(3)
	faddl	-56(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fldl	56(%eax)
	fxch	%st(5)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-40(%eax)
	fldl	40(%eax)
	fxch	%st(6)
	fstpl	64(%eax)
	fxch	%st(1)
	fstl	-24(%eax)
	fxch	%st(5)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(5)
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(4), %st
	fxch	%st(4)
	fstpl	-8(%eax)
	fxch	%st(3)
	fstl	8(%eax)
	fldl	72(%eax)
	fldl	88(%eax)
	fldl	120(%eax)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(5), %st
	fxch	%st(5)
	faddp	%st, %st(3)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	24(%eax)
	fldl	104(%eax)
	fxch	%st(2)
	fstl	40(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(4)
	faddp	%st, %st(2)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	56(%eax)
	fxch	%st(1)
	fstpl	120(%eax)
	fxch	%st(2)
	fstpl	104(%eax)
	fstpl	88(%eax)
	fstpl	72(%eax)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fldl	-48(%edx)
	fldl	(%edx)
	fldl	16(%edx)
	fldl	48(%edx)
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fldl	64(%edx)
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(2)
	fadd	%st(6), %st
	fxch	%st(4)
	faddl	-64(%edx)
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(6), %st
	fxch	%st(6)
	fstpl	-48(%edx)
	fldl	32(%edx)
	fxch	%st(6)
	fstl	-32(%edx)
	fxch	%st(6)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(6)
	fxch	%st(2)
	fadd	%st(3), %st
	fldl	112(%edx)
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(6)
	fstl	-16(%edx)
	faddp	%st, %st(2)
	fadd	%st, %st(4)
	fldl	80(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(6), %st
	fxch	%st(2)
	fstl	(%edx)
	faddp	%st, %st(6)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(5), %st
	fxch	%st(5)
	fstpl	16(%edx)
	fldl	96(%edx)
	fxch	%st(5)
	fstl	32(%edx)
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(5)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fstpl	48(%edx)
	fxch	%st(2)
	fstpl	112(%edx)
	fxch	%st(1)
	fstpl	96(%edx)
	fstpl	80(%edx)
	fstpl	64(%edx)
	fldl	-8(%edx)
	fldl	-24(%edx)
	fldl	-40(%edx)
	fldl	8(%edx)
	fldl	24(%edx)
	fldl	56(%edx)
	fxch	%st(2)
	fadd	%st(5), %st
	fxch	%st(5)
	fadd	%st(4), %st
	fldl	72(%edx)
	fxch	%st(5)
	fadd	%st(4), %st
	fxch	%st(2)
	cmpl	%ecx, -1392(%ebp)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(1), %st
	fxch	%st(4)
	faddl	-56(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-40(%edx)
	fldl	40(%edx)
	fxch	%st(1)
	fstl	-24(%edx)
	fxch	%st(1)
	fadd	%st(6), %st
	fxch	%st(6)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(2)
	fadd	%st(5), %st
	fldl	120(%edx)
	fxch	%st(6)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(3)
	fstl	-8(%edx)
	faddp	%st, %st(2)
	fadd	%st, %st(4)
	fldl	88(%edx)
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(2)
	fstl	8(%edx)
	faddp	%st, %st(3)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	24(%edx)
	fldl	104(%edx)
	fxch	%st(2)
	fstl	40(%edx)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(2)
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	56(%edx)
	fxch	%st(2)
	fstpl	120(%edx)
	fstpl	104(%edx)
	fxch	%st(1)
	fstpl	88(%edx)
	fstpl	72(%edx)
L634:
	jge	L643
	movl	-1404(%ebp), %edx
	subl	$3, %edx
	cmpl	%ecx, %edx
	jl	L405
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$4, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	-32(%eax)
	fldl	-48(%eax)
	fldl	(%eax)
	fldl	16(%eax)
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(2)
	faddl	-64(%eax)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fstpl	-48(%eax)
	fldl	32(%eax)
	fxch	%st(3)
	fstl	-32(%eax)
	fxch	%st(3)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(3)
	fldl	48(%eax)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(3)
	fstl	-16(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(2)
	fxch	%st(3)
	fstpl	48(%eax)
	fxch	%st(1)
	fstpl	32(%eax)
	fldl	-8(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fldl	-24(%eax)
	fldl	-40(%eax)
	fldl	8(%eax)
	fxch	%st(4)
	fstpl	16(%eax)
	fldl	24(%eax)
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(1)
	faddl	-56(%eax)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-40(%eax)
	fldl	40(%eax)
	fxch	%st(2)
	fstl	-24(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	faddp	%st, %st(2)
	fldl	56(%eax)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(2)
	fstl	-8(%eax)
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	faddp	%st, %st(1)
	fxch	%st(2)
	fstpl	56(%eax)
	fstpl	40(%eax)
	fxch	%st(1)
	fstpl	24(%eax)
	fstpl	8(%eax)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fldl	-48(%edx)
	fldl	(%edx)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-64(%edx)
	fstl	-48(%edx)
	fldl	16(%edx)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(1)
	fldl	32(%edx)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(1)
	fstl	-32(%edx)
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(3)
	fldl	48(%edx)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(3)
	fstl	-16(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(3)
	fstpl	48(%edx)
	fldl	-8(%edx)
	fldl	-24(%edx)
	fxch	%st(3)
	fstpl	32(%edx)
	fxch	%st(3)
	fstpl	(%edx)
	fldl	-40(%edx)
	fldl	8(%edx)
	fxch	%st(2)
	fstpl	16(%edx)
	fldl	24(%edx)
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(1)
	faddl	-56(%edx)
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fstpl	-40(%edx)
	fldl	40(%edx)
	fxch	%st(3)
	fstl	-24(%edx)
	fxch	%st(3)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(3)
	fldl	56(%edx)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(3)
	fstl	-8(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(3)
	fstpl	56(%edx)
	fxch	%st(1)
	fstpl	40(%edx)
	fstpl	24(%edx)
	fstpl	8(%edx)
L405:
	movl	-1404(%ebp), %edx
	decl	%edx
	cmpl	%ecx, %edx
	jl	L406
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$2, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	-32(%eax)
	fldl	-48(%eax)
	fldl	(%eax)
	fldl	16(%eax)
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(2)
	faddl	-64(%eax)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(4)
	fstpl	16(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-48(%eax)
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	-40(%eax)
	fxch	%st(1)
	fstpl	-32(%eax)
	fldl	-8(%eax)
	fldl	-24(%eax)
	fldl	8(%eax)
	fldl	24(%eax)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(4), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(4)
	faddl	-56(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(3)
	fstpl	24(%eax)
	fxch	%st(3)
	fstpl	8(%eax)
	fadd	%st, %st(2)
	fstpl	-40(%eax)
	fstpl	-8(%eax)
	fstpl	-24(%eax)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fldl	-48(%edx)
	fldl	(%edx)
	fldl	16(%edx)
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(1)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(3), %st
	fxch	%st(2)
	faddl	-64(%edx)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(4)
	fstpl	16(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-48(%edx)
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	-40(%edx)
	fldl	-8(%edx)
	fxch	%st(2)
	fstpl	-32(%edx)
	fldl	8(%edx)
	fldl	-24(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	-56(%edx)
	fstl	-40(%edx)
	fldl	24(%edx)
	fadd	%st(4), %st
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fxch	%st(3)
	fstpl	24(%edx)
	fstpl	8(%edx)
	fstpl	-8(%edx)
	fstpl	-24(%edx)
L406:
	cmpl	-1404(%ebp), %ecx
	jg	L407
	movl	%ecx, %edx
	sall	$4, %edx
	incl	%ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	(%eax)
	fadd	%st(1), %st
	fstpl	(%eax)
	fldl	-32(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	-48(%eax)
	fadd	%st, %st(1)
	faddl	-64(%eax)
	fxch	%st(1)
	fstpl	-32(%eax)
	fldl	8(%eax)
	fxch	%st(1)
	fstpl	-48(%eax)
	fldl	-8(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	8(%eax)
	fldl	-24(%eax)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-8(%eax)
	fldl	-40(%eax)
	fadd	%st, %st(1)
	faddl	-56(%eax)
	fxch	%st(1)
	fstpl	-24(%eax)
	fstpl	-40(%eax)
	fldl	-16(%edx)
	fldl	(%edx)
	fadd	%st(1), %st
	fstpl	(%edx)
	fldl	-32(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	-48(%edx)
	fadd	%st, %st(1)
	faddl	-64(%edx)
	fxch	%st(1)
	fstpl	-32(%edx)
	fldl	8(%edx)
	fxch	%st(1)
	fstpl	-48(%edx)
	fldl	-8(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	8(%edx)
	fldl	-24(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-8(%edx)
	fldl	-40(%edx)
	fadd	%st, %st(1)
	fxch	%st(1)
	fstpl	-24(%edx)
	faddl	-56(%edx)
	fstpl	-40(%edx)
L407:
	leal	-1(%ecx), %edx
	subl	$4, %edi
	sall	$4, %edx
	leal	(%edx,%ebx), %ecx
	addl	%esi, %edx
	cmpl	$3, %edi
	fldl	-16(%ecx)
	fldl	-32(%ecx)
	fldl	(%ecx)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-48(%ecx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-32(%ecx)
	fldl	-24(%ecx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-16(%ecx)
	fldl	-8(%ecx)
	fxch	%st(2)
	fstpl	(%ecx)
	fldl	8(%ecx)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-40(%ecx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-24(%ecx)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-8(%ecx)
	fstpl	8(%ecx)
	fldl	-16(%edx)
	fldl	-32(%edx)
	fldl	(%edx)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-48(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-32(%edx)
	fldl	-24(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-16(%edx)
	fldl	-8(%edx)
	fxch	%st(2)
	fstpl	(%edx)
	fldl	8(%edx)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-40(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-24(%edx)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-8(%edx)
	fstpl	8(%edx)
	jg	L408
L595:
	cmpl	$1, %edi
	jle	L409
	movl	-1404(%ebp), %eax
	movl	%edi, %ecx
	subl	$7, %eax
	movl	%eax, -1396(%ebp)
	cmpl	%edi, %eax
	jmp	L635
L644:
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$8, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	(%eax)
	fldl	16(%eax)
	fldl	80(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	faddl	-32(%eax)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	32(%eax)
	fxch	%st(1)
	fstl	(%eax)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	48(%eax)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%eax)
	fldl	64(%eax)
	fxch	%st(1)
	fstl	32(%eax)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	112(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	48(%eax)
	fldl	96(%eax)
	fxch	%st(1)
	fstl	64(%eax)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	24(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	80(%eax)
	fldl	88(%eax)
	fxch	%st(3)
	fstpl	112(%eax)
	fstpl	96(%eax)
	fldl	-8(%eax)
	fldl	8(%eax)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-24(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-8(%eax)
	fldl	40(%eax)
	fxch	%st(2)
	fstl	8(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(2)
	fldl	56(%eax)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	24(%eax)
	fldl	72(%eax)
	fxch	%st(2)
	fstl	40(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(2)
	fadd	%st, %st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	56(%eax)
	fldl	104(%eax)
	fxch	%st(1)
	fstl	72(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(1)
	fstl	88(%eax)
	fldl	120(%eax)
	fadd	%st(2), %st
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstpl	120(%eax)
	fstpl	104(%eax)
	fldl	-16(%edx)
	fldl	(%edx)
	fldl	16(%edx)
	fldl	80(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	faddl	-32(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	32(%edx)
	fxch	%st(1)
	fstl	(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	48(%edx)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%edx)
	fldl	64(%edx)
	fxch	%st(1)
	fstl	32(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	112(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	48(%edx)
	fldl	96(%edx)
	fxch	%st(1)
	fstl	64(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	24(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	80(%edx)
	fxch	%st(2)
	fstpl	112(%edx)
	fxch	%st(1)
	fstpl	96(%edx)
	fldl	-8(%edx)
	fldl	8(%edx)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-24(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	-8(%edx)
	fldl	40(%edx)
	fxch	%st(2)
	fstl	8(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(2)
	fldl	56(%edx)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	24(%edx)
	fldl	72(%edx)
	fxch	%st(2)
	fstl	40(%edx)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddp	%st, %st(2)
	fxch	%st(1)
	fstl	56(%edx)
	fldl	88(%edx)
	fldl	120(%edx)
	fxch	%st(1)
	cmpl	%ecx, -1396(%ebp)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(2)
	fldl	104(%edx)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	72(%edx)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	88(%edx)
	fxch	%st(1)
	fstpl	120(%edx)
	fstpl	104(%edx)
L635:
	jge	L644
	movl	-1404(%ebp), %edx
	subl	$3, %edx
	cmpl	%ecx, %edx
	jl	L414
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$4, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	(%eax)
	fldl	16(%eax)
	fldl	32(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	faddl	-32(%eax)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-16(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fldl	48(%eax)
	fxch	%st(2)
	fstpl	(%eax)
	fldl	8(%eax)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%eax)
	fldl	24(%eax)
	fxch	%st(3)
	fstpl	48(%eax)
	fldl	-8(%eax)
	fxch	%st(1)
	fstpl	32(%eax)
	fldl	40(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-24(%eax)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fstpl	-8(%eax)
	fadd	%st, %st(1)
	fadd	%st(2), %st
	fldl	56(%eax)
	fxch	%st(3)
	fstpl	8(%eax)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	24(%eax)
	fstpl	56(%eax)
	fstpl	40(%eax)
	fldl	-16(%edx)
	fldl	(%edx)
	fldl	16(%edx)
	fldl	32(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	faddl	-32(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-16(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fldl	48(%edx)
	fxch	%st(2)
	fstpl	(%edx)
	fldl	8(%edx)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	16(%edx)
	fxch	%st(2)
	fstpl	48(%edx)
	fldl	-8(%edx)
	fxch	%st(2)
	fstpl	32(%edx)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-24(%edx)
	fstl	-8(%edx)
	fldl	24(%edx)
	fldl	40(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	faddp	%st, %st(2)
	fldl	56(%edx)
	fxch	%st(1)
	fadd	%st(3), %st
	fxch	%st(3)
	fadd	%st(2), %st
	fxch	%st(2)
	fstpl	8(%edx)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	24(%edx)
	fxch	%st(1)
	fstpl	56(%edx)
	fstpl	40(%edx)
L414:
	movl	-1404(%ebp), %edx
	decl	%edx
	cmpl	%ecx, %edx
	jl	L415
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$2, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	(%eax)
	fldl	16(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	-32(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	-8(%eax)
	fxch	%st(2)
	fstpl	16(%eax)
	fstpl	(%eax)
	fldl	8(%eax)
	fldl	24(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	-24(%eax)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-8(%eax)
	fxch	%st(1)
	fstpl	24(%eax)
	fstpl	8(%eax)
	fldl	-16(%edx)
	fldl	(%edx)
	fldl	16(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	-32(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	-8(%edx)
	fxch	%st(2)
	fstpl	16(%edx)
	fstpl	(%edx)
	fldl	8(%edx)
	fldl	24(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	faddl	-24(%edx)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstpl	-8(%edx)
	fxch	%st(1)
	fstpl	24(%edx)
	fstpl	8(%edx)
L415:
	cmpl	-1404(%ebp), %ecx
	jg	L416
	movl	%ecx, %edx
	sall	$4, %edx
	incl	%ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	fldl	(%eax)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-32(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fldl	8(%eax)
	fxch	%st(1)
	fstpl	-16(%eax)
	fldl	-8(%eax)
	fadd	%st, %st(1)
	faddl	-24(%eax)
	fxch	%st(1)
	fstpl	8(%eax)
	fstpl	-8(%eax)
	fldl	-16(%edx)
	fldl	(%edx)
	fadd	%st(1), %st
	fxch	%st(1)
	faddl	-32(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fldl	8(%edx)
	fxch	%st(1)
	fstpl	-16(%edx)
	fldl	-8(%edx)
	fadd	%st, %st(1)
	faddl	-24(%edx)
	fxch	%st(1)
	fstpl	8(%edx)
	fstpl	-8(%edx)
L416:
	leal	-1(%ecx), %eax
	subl	$2, %edi
	sall	$4, %eax
	leal	(%eax,%ebx), %ecx
	addl	%esi, %eax
	fldl	-16(%ecx)
	faddl	(%ecx)
	fstpl	(%ecx)
	fldl	-8(%ecx)
	faddl	8(%ecx)
	fstpl	8(%ecx)
	fldl	-16(%eax)
	faddl	(%eax)
	fstpl	(%eax)
	fldl	-8(%eax)
	faddl	8(%eax)
	fstpl	8(%eax)
L409:
	testl	%edi, %edi
	jle	L417
	movl	%edi, %ecx
	movl	-1404(%ebp), %edi
	subl	$7, %edi
	jmp	L636
	.p2align 4,,7
L645:
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$8, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	faddl	(%eax)
	fstl	(%eax)
	faddl	16(%eax)
	fstl	16(%eax)
	faddl	32(%eax)
	fstl	32(%eax)
	faddl	48(%eax)
	fstl	48(%eax)
	faddl	64(%eax)
	fstl	64(%eax)
	faddl	80(%eax)
	fstl	80(%eax)
	faddl	96(%eax)
	fstl	96(%eax)
	faddl	112(%eax)
	fstpl	112(%eax)
	fldl	-8(%eax)
	faddl	8(%eax)
	fstl	8(%eax)
	faddl	24(%eax)
	fstl	24(%eax)
	faddl	40(%eax)
	fstl	40(%eax)
	faddl	56(%eax)
	fstl	56(%eax)
	faddl	72(%eax)
	fstl	72(%eax)
	faddl	88(%eax)
	fstl	88(%eax)
	faddl	104(%eax)
	fstl	104(%eax)
	faddl	120(%eax)
	fstpl	120(%eax)
	fldl	-16(%edx)
	faddl	(%edx)
	fstl	(%edx)
	faddl	16(%edx)
	fstl	16(%edx)
	faddl	32(%edx)
	fstl	32(%edx)
	faddl	48(%edx)
	fstl	48(%edx)
	faddl	64(%edx)
	fstl	64(%edx)
	faddl	80(%edx)
	fstl	80(%edx)
	faddl	96(%edx)
	fstl	96(%edx)
	faddl	112(%edx)
	fstpl	112(%edx)
	fldl	-8(%edx)
	faddl	8(%edx)
	fstl	8(%edx)
	faddl	24(%edx)
	fstl	24(%edx)
	faddl	40(%edx)
	fstl	40(%edx)
	faddl	56(%edx)
	fstl	56(%edx)
	faddl	72(%edx)
	fstl	72(%edx)
	faddl	88(%edx)
	fstl	88(%edx)
	faddl	104(%edx)
	fstl	104(%edx)
	faddl	120(%edx)
	fstpl	120(%edx)
L636:
	cmpl	%ecx, %edi
	jge	L645
	movl	-1404(%ebp), %edi
	subl	$3, %edi
	cmpl	%ecx, %edi
	jl	L422
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$4, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	faddl	(%eax)
	fstl	(%eax)
	faddl	16(%eax)
	fstl	16(%eax)
	faddl	32(%eax)
	fstl	32(%eax)
	faddl	48(%eax)
	fstpl	48(%eax)
	fldl	-8(%eax)
	faddl	8(%eax)
	fstl	8(%eax)
	faddl	24(%eax)
	fstl	24(%eax)
	faddl	40(%eax)
	fstl	40(%eax)
	faddl	56(%eax)
	fstpl	56(%eax)
	fldl	-16(%edx)
	faddl	(%edx)
	fstl	(%edx)
	faddl	16(%edx)
	fstl	16(%edx)
	faddl	32(%edx)
	fstl	32(%edx)
	faddl	48(%edx)
	fstpl	48(%edx)
	fldl	-8(%edx)
	faddl	8(%edx)
	fstl	8(%edx)
	faddl	24(%edx)
	fstl	24(%edx)
	faddl	40(%edx)
	fstl	40(%edx)
	faddl	56(%edx)
	fstpl	56(%edx)
L422:
	movl	-1404(%ebp), %edi
	decl	%edi
	cmpl	%ecx, %edi
	jl	L423
	movl	%ecx, %edx
	sall	$4, %edx
	addl	$2, %ecx
	leal	(%edx,%ebx), %eax
	addl	%esi, %edx
	fldl	-16(%eax)
	faddl	(%eax)
	fstl	(%eax)
	faddl	16(%eax)
	fstpl	16(%eax)
	fldl	-8(%eax)
	faddl	8(%eax)
	fstl	8(%eax)
	faddl	24(%eax)
	fstpl	24(%eax)
	fldl	-16(%edx)
	faddl	(%edx)
	fstl	(%edx)
	faddl	16(%edx)
	fstpl	16(%edx)
	fldl	-8(%edx)
	faddl	8(%edx)
	fstl	8(%edx)
	faddl	24(%edx)
	fstpl	24(%edx)
L423:
	cmpl	-1404(%ebp), %ecx
	jg	L417
	movl	%ecx, %edi
	sall	$4, %edi
	leal	(%edi,%ebx), %ecx
	fldl	-16(%ecx)
	addl	%esi, %edi
	faddl	(%ecx)
	fstpl	(%ecx)
	fldl	-8(%ecx)
	faddl	8(%ecx)
	fstpl	8(%ecx)
	fldl	-16(%edi)
	faddl	(%edi)
	fstpl	(%edi)
	fldl	-8(%edi)
	faddl	8(%edi)
	fstpl	8(%edi)
L417:
	movl	_potshift, %edx
	testl	%edx, %edx
	jne	L425
	sall	$4, -1404(%ebp)
	movl	-1400(%ebp), %ecx
	movl	56(%ebp), %edi
	fldl	_I
	movl	-1404(%ebp), %edx
	movl	44(%ebp), %eax
	fld	%st(0)
	addl	$16, %edx
	imull	%ecx, %edx
	fldl	(%eax)
	movl	52(%ebp), %ecx
	fldl	8(%eax)
	addl	%edx, %edi
	addl	%ecx, %edx
	fldl	(%edx)
	fxch	%st(2)
	movl	%edi, -1312(%ebp)
	movl	-1300(%ebp), %edi
	fstpl	-1320(%ebp)
	fxch	%st(1)
	fstpl	-1336(%ebp)
	fldl	_I+8
	fldl	8(%edx)
	fxch	%st(4)
	fstl	-144(%ebp)
	fxch	%st(2)
	fstpl	-1328(%ebp)
	fstl	-136(%ebp)
	fldl	24(%edi)
	fxch	%st(4)
	fstpl	-1344(%ebp)
	fldz
	fxch	%st(3)
	fmul	%st(4), %st
	fxch	%st(3)
	fmul	%st(1), %st
	fxch	%st(4)
	fmulp	%st, %st(1)
	fxch	%st(2)
	fsubp	%st, %st(3)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(2)
	fstl	-160(%ebp)
	fstl	-144(%ebp)
	fstl	-128(%ebp)
	fxch	%st(2)
	faddp	%st, %st(1)
	fld	%st(0)
	fstl	-136(%ebp)
	fxch	%st(1)
	fabs
	fxch	%st(1)
	fstl	-120(%ebp)
	fstl	-152(%ebp)
	fxch	%st(2)
	faddl	16(%edi)
	fxch	%st(2)
	leal	-112(%ebp), %edi
	fstl	-104(%ebp)
	fld	%st(2)
	fabs
	fxch	%st(3)
	fstl	-160(%ebp)
	fstl	-112(%ebp)
	fxch	%st(2)
	fucom	%st(3)
	fnstsw	%ax
	fxch	%st(3)
	fstpl	-1272(%ebp)
	fxch	%st(2)
	sahf
	fstpl	-1280(%ebp)
	jbe	L438
	leal	-1280(%ebp), %eax
L440:
	fldl	(%eax)
	fldz
	fld	%st(1)
	fucom	%st(1)
	fnstsw	%ax
	fstp	%st(1)
	sahf
	jp	L651
	je	L652
L651:
	fstp	%st(1)
	fdivr	%st, %st(1)
	fdivr	%st, %st(2)
	fxch	%st(1)
	fmul	%st(0), %st
	fxch	%st(2)
	fmul	%st(0), %st
	faddp	%st, %st(2)
	fxch	%st(1)
	fsqrt
	fmulp	%st, %st(1)
L443:
	fstpl	(%esp)
	call	_log
	fstpl	-1352(%ebp)
	fldl	(%edi)
	fstpl	8(%esp)
	fldl	8(%edi)
	fstpl	(%esp)
	call	_atan2
	fldl	-1352(%ebp)
	fld	%st(1)
	fldl	_I
	fld	%st(2)
	fldl	_I+8
	fld	%st(2)
	fstl	-240(%ebp)
	fxch	%st(5)
	fstl	-176(%ebp)
	fxch	%st(6)
	fstl	-168(%ebp)
	fxch	%st(6)
	fstpl	-96(%ebp)
	fldz
	fmul	%st(5), %st
	fxch	%st(6)
	fstpl	-88(%ebp)
	fldz
	fxch	%st(1)
	fstpl	-232(%ebp)
	fldl	8(%ebx)
	fxch	%st(1)
	fmull	_I+8
	fxch	%st(3)
	fmul	%st(1), %st
	fldl	-1336(%ebp)
	fxch	%st(2)
	fmull	_I+8
	fxch	%st(1)
	fsub	%st(4), %st
	fxch	%st(2)
	fmul	%st(3), %st
	fxch	%st(4)
	fstpl	-1360(%ebp)
	fadd	%st(6), %st
	fldl	-1344(%ebp)
	fxch	%st(7)
	fstpl	-1368(%ebp)
	fldl	-1336(%ebp)
	fxch	%st(7)
	fmul	%st(5), %st
	fxch	%st(7)
	fmul	%st(5), %st
	fxch	%st(1)
	fstl	-232(%ebp)
	fxch	%st(4)
	fsubp	%st, %st(7)
	fldl	-1344(%ebp)
	fxch	%st(4)
	fstl	-216(%ebp)
	fxch	%st(2)
	fstl	-256(%ebp)
	fxch	%st(4)
	fmul	%st(3), %st
	fxch	%st(2)
	fstl	-248(%ebp)
	fxch	%st(4)
	fstl	-240(%ebp)
	fstl	-224(%ebp)
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(1)
	faddl	(%ebx)
	fxch	%st(3)
	fstl	-200(%ebp)
	fxch	%st(5)
	fmul	%st(1), %st
	fxch	%st(1)
	fmull	_I+8
	fxch	%st(1)
	fsubl	-1360(%ebp)
	fxch	%st(3)
	fstl	-256(%ebp)
	fstl	-208(%ebp)
	fxch	%st(1)
	faddl	-1368(%ebp)
	fxch	%st(6)
	fadd	%st(3), %st
	fxch	%st(3)
	fstl	-304(%ebp)
	fxch	%st(1)
	fsub	%st(3), %st
	fxch	%st(5)
	fsub	%st(6), %st
	fxch	%st(6)
	fstl	-296(%ebp)
	fxch	%st(1)
	fstpl	-288(%ebp)
	movl	48(%ebp), %ecx
	fstl	-280(%ebp)
	fstl	-312(%ebp)
	fxch	%st(2)
	fstl	-320(%ebp)
	fstpl	-272(%ebp)
	fxch	%st(1)
	fstpl	-264(%ebp)
	fldz
	fxch	%st(3)
	fstl	-336(%ebp)
	fxch	%st(4)
	fstl	-328(%ebp)
	fxch	%st(4)
	fstl	-192(%ebp)
	fxch	%st(4)
	fstl	-184(%ebp)
	fxch	%st(4)
	faddl	(%ecx)
	fxch	%st(4)
	faddl	8(%ecx)
	fxch	%st(2)
	fsubl	LC18
	fxch	%st(4)
	fstpl	(%ecx)
	fxch	%st(1)
	fstpl	8(%ecx)
	fldz
	fldl	_I
	fldl	_I+8
	fxch	%st(2)
	fmul	%st(1), %st
	fxch	%st(4)
	fmul	%st(2), %st
	fxch	%st(1)
	fstl	-400(%ebp)
	fld	%st(0)
	fxch	%st(3)
	fstl	-392(%ebp)
	fldl	8(%esi)
	fxch	%st(3)
	fstl	-1376(%ebp)
	fxch	%st(6)
	fstl	-1384(%ebp)
	fxch	%st(4)
	fmul	%st(3), %st
	fxch	%st(3)
	fmul	%st(1), %st
	fxch	%st(3)
	fsubp	%st, %st(6)
	fldl	-1328(%ebp)
	fxch	%st(3)
	faddp	%st, %st(4)
	fldl	-1320(%ebp)
	fxch	%st(3)
	fmul	%st(7), %st
	fxch	%st(3)
	fmul	%st(5), %st
	fxch	%st(7)
	fmull	-1320(%ebp)
	fxch	%st(5)
	fmull	-1328(%ebp)
	fxch	%st(7)
	fsubp	%st, %st(3)
	fxch	%st(3)
	fstl	-392(%ebp)
	fstl	-376(%ebp)
	fxch	%st(5)
	fstl	-416(%ebp)
	fxch	%st(6)
	faddp	%st, %st(4)
	fxch	%st(4)
	fstl	-408(%ebp)
	fxch	%st(5)
	fstl	-400(%ebp)
	fxch	%st(4)
	fmul	%st(3), %st
	fxch	%st(3)
	fmulp	%st, %st(2)
	fxch	%st(2)
	fsubl	-1376(%ebp)
	fxch	%st(3)
	fstl	-384(%ebp)
	faddl	(%esi)
	fxch	%st(1)
	faddl	-1384(%ebp)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(4)
	fstl	-360(%ebp)
	fsub	%st(2), %st
	fxch	%st(1)
	fstl	-416(%ebp)
	fstl	-368(%ebp)
	fsub	%st(4), %st
	fxch	%st(3)
	fstl	-464(%ebp)
	fxch	%st(2)
	movl	-1312(%ebp), %edi
	fstl	-456(%ebp)
	fxch	%st(2)
	movl	_pshift, %edx
	fstpl	-448(%ebp)
	fxch	%st(1)
	fstl	-440(%ebp)
	fstl	-472(%ebp)
	fxch	%st(3)
	fstl	-480(%ebp)
	fstpl	-432(%ebp)
	fxch	%st(2)
	fstpl	-424(%ebp)
	fstl	-496(%ebp)
	fxch	%st(1)
	fstl	-488(%ebp)
	fxch	%st(1)
	fstl	-352(%ebp)
	fxch	%st(1)
	fstl	-344(%ebp)
	fxch	%st(1)
	faddl	(%edi)
	fxch	%st(1)
	faddl	8(%edi)
	fxch	%st(1)
	movl	%edx, -1404(%ebp)
	fstpl	(%edi)
	fstpl	8(%edi)
	movl	$1, %edi
	cmpl	%edx, %edi
	jmp	L637
L646:
	pushl	%edi
	fldl	-1336(%ebp)
	movl	%edi, %eax
	movl	-1300(%ebp), %ecx
	fildl	(%esp)
	sall	$4, %eax
	fldl	-1344(%ebp)
	fldl	(%eax,%ecx)
	addl	%eax, %ecx
	fldl	8(%ecx)
	fxch	%st(4)
	movl	%ecx, -1416(%ebp)
	movl	48(%ebp), %ecx
	fdiv	%st(3), %st
	fxch	%st(2)
	leal	(%eax,%ecx), %edx
	fdiv	%st(3), %st
	fxch	%st(2)
	fsubrl	(%eax,%ebx)
	fxch	%st(2)
	fsubrl	8(%eax,%ebx)
	fld	%st(2)
	fld	%st(1)
	fxch	%st(4)
	fmul	%st(6), %st
	fxch	%st(2)
	fmul	%st(3), %st
	fxch	%st(4)
	fmulp	%st, %st(6)
	fmulp	%st, %st(2)
	faddp	%st, %st(2)
	fldl	_I
	fldl	_I+8
	fxch	%st(2)
	fsubp	%st, %st(5)
	fld	%st(0)
	fldz
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(4)
	fmulp	%st, %st(3)
	fsubrp	%st, %st(3)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(4)
	fadd	%st(2), %st
	fxch	%st(2)
	fstl	-544(%ebp)
	fstpl	-528(%ebp)
	faddp	%st, %st(3)
	fldl	-1320(%ebp)
	fxch	%st(1)
	fstl	-560(%ebp)
	fstl	-512(%ebp)
	fxch	%st(1)
	fdiv	%st(2), %st
	fxch	%st(2)
	fdivrl	-1328(%ebp)
	fxch	%st(3)
	fstl	-536(%ebp)
	fstl	-520(%ebp)
	fstl	-552(%ebp)
	fstl	-504(%ebp)
	fxch	%st(1)
	faddl	(%edx)
	fxch	%st(1)
	faddl	8(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fstpl	8(%edx)
	movl	-1416(%ebp), %edx
	fsubrl	(%eax,%esi)
	fxch	%st(1)
	fsubrl	8(%eax,%esi)
	fldl	(%edx)
	fldl	8(%edx)
	fld	%st(3)
	movl	-1312(%ebp), %edx
	fld	%st(3)
	fxch	%st(5)
	fmul	%st(2), %st
	fxch	%st(4)
	addl	%edx, %eax
	fmul	%st(3), %st
	fxch	%st(5)
	fmulp	%st, %st(2)
	fmulp	%st, %st(2)
	fldl	_I
	fxch	%st(3)
	faddp	%st, %st(4)
	fldl	_I+8
	fxch	%st(2)
	fsubp	%st, %st(1)
	fld	%st(2)
	fldz
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(5)
	fmulp	%st, %st(3)
	fsubrp	%st, %st(4)
	fldz
	fmulp	%st, %st(3)
	fadd	%st(3), %st
	fxch	%st(3)
	fstl	-608(%ebp)
	fstpl	-592(%ebp)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-624(%ebp)
	fstl	-576(%ebp)
	fldl	-1336(%ebp)
	fxch	%st(2)
	fstl	-600(%ebp)
	fstl	-584(%ebp)
	fstl	-616(%ebp)
	fstl	-568(%ebp)
	fxch	%st(1)
	fsubrl	(%eax)
	fxch	%st(1)
	fsubrl	8(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fstpl	8(%eax)
	leal	1(%edi), %eax
	movl	%eax, %edx
	movl	%eax, (%esp)
	fldl	-1344(%ebp)
	addl	$2, %edi
	movl	-1300(%ebp), %eax
	fildl	(%esp)
	addl	$4, %esp
	sall	$4, %edx
	fldl	(%edx,%eax)
	addl	%edx, %eax
	fldl	8(%eax)
	fxch	%st(4)
	movl	%eax, -1416(%ebp)
	leal	(%edx,%ecx), %eax
	fdiv	%st(2), %st
	fsubrl	(%edx,%ebx)
	fxch	%st(3)
	fdiv	%st(2), %st
	fld	%st(3)
	fxch	%st(1)
	fsubrl	8(%edx,%ebx)
	fxch	%st(4)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(2), %st
	fld	%st(4)
	fmulp	%st, %st(3)
	fxch	%st(4)
	fmulp	%st, %st(5)
	fldl	_I+8
	fxch	%st(1)
	faddp	%st, %st(2)
	fldl	_I
	fxch	%st(4)
	fsubp	%st, %st(5)
	fldz
	fld	%st(4)
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(2), %st
	fxch	%st(3)
	fmulp	%st, %st(2)
	fsubp	%st, %st(2)
	fldz
	fmulp	%st, %st(4)
	fxch	%st(4)
	fadd	%st(1), %st
	fxch	%st(1)
	fstl	-672(%ebp)
	fstpl	-656(%ebp)
	fxch	%st(2)
	faddp	%st, %st(3)
	fxch	%st(1)
	fstl	-688(%ebp)
	fstl	-640(%ebp)
	fxch	%st(2)
	fstl	-664(%ebp)
	fstl	-648(%ebp)
	fstl	-680(%ebp)
	fstl	-632(%ebp)
	fxch	%st(2)
	faddl	(%eax)
	fxch	%st(2)
	faddl	8(%eax)
	fxch	%st(2)
	fstpl	(%eax)
	fxch	%st(1)
	fstpl	8(%eax)
	movl	-1416(%ebp), %ecx
	movl	-1312(%ebp), %eax
	fldl	-1320(%ebp)
	fldl	(%ecx)
	fldl	8(%ecx)
	fxch	%st(2)
	fdiv	%st(3), %st
	fxch	%st(3)
	fdivrl	-1328(%ebp)
	fxch	%st(3)
	fsubrl	(%edx,%esi)
	fxch	%st(3)
	fsubrl	8(%edx,%esi)
	addl	%eax, %edx
	cmpl	-1404(%ebp), %edi
	fld	%st(3)
	fld	%st(1)
	fxch	%st(5)
	fmul	%st(4), %st
	fxch	%st(2)
	fmul	%st(3), %st
	fxch	%st(5)
	fmulp	%st, %st(4)
	fldl	_I
	fxch	%st(1)
	fmulp	%st, %st(3)
	fxch	%st(1)
	faddp	%st, %st(4)
	fldl	_I+8
	fld	%st(1)
	fxch	%st(3)
	fsubp	%st, %st(4)
	fldz
	fxch	%st(3)
	fmul	%st(5), %st
	fxch	%st(3)
	fmul	%st(1), %st
	fxch	%st(5)
	fmulp	%st, %st(1)
	fxch	%st(2)
	fsubp	%st, %st(4)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fstl	-736(%ebp)
	fstpl	-720(%ebp)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-752(%ebp)
	fstl	-704(%ebp)
	fxch	%st(1)
	fstl	-728(%ebp)
	fstl	-712(%ebp)
	fstl	-744(%ebp)
	fstl	-696(%ebp)
	fxch	%st(1)
	faddl	(%edx)
	fxch	%st(1)
	faddl	8(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fstpl	8(%edx)
L637:
	jl	L646
	testb	$1, -1404(%ebp)
	je	L606
	fildl	-1404(%ebp)
	movl	-1404(%ebp), %edi
	movl	-1300(%ebp), %ecx
	fldl	-1336(%ebp)
	sall	$4, %edi
	leal	(%edi,%ecx), %eax
	fldl	-1344(%ebp)
	fldl	(%edi,%ecx)
	fxch	%st(2)
	fdiv	%st(3), %st
	fxch	%st(1)
	fdiv	%st(3), %st
	fldl	8(%eax)
	fxch	%st(2)
	fsubrl	(%edi,%ebx)
	fxch	%st(1)
	fsubrl	8(%edi,%ebx)
	movl	48(%ebp), %ebx
	fld	%st(1)
	fld	%st(1)
	fxch	%st(3)
	fmul	%st(4), %st
	fxch	%st(2)
	leal	(%edi,%ebx), %edx
	fmul	%st(5), %st
	fxch	%st(3)
	fmulp	%st, %st(4)
	fmulp	%st, %st(4)
	fldl	_I+8
	fxch	%st(1)
	faddp	%st, %st(2)
	fldl	_I
	fxch	%st(4)
	fsubp	%st, %st(3)
	fldz
	fld	%st(4)
	fxch	%st(1)
	fmul	%st(2), %st
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(3)
	fmulp	%st, %st(2)
	fsubrp	%st, %st(2)
	fldz
	fmulp	%st, %st(4)
	fxch	%st(2)
	fadd	%st(1), %st
	fxch	%st(1)
	fstl	-800(%ebp)
	fstpl	-784(%ebp)
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-816(%ebp)
	fstl	-768(%ebp)
	fxch	%st(1)
	fstl	-792(%ebp)
	fstl	-776(%ebp)
	fstl	-808(%ebp)
	fstl	-760(%ebp)
	faddl	8(%edx)
	fxch	%st(1)
	faddl	(%edx)
	fxch	%st(1)
	fstpl	8(%edx)
	fldl	-1320(%ebp)
	fxch	%st(1)
	fstpl	(%edx)
	fldl	(%eax)
	fxch	%st(1)
	fdiv	%st(2), %st
	fxch	%st(2)
	fdivrl	-1328(%ebp)
	fxch	%st(2)
	fsubrl	(%edi,%esi)
	fld	%st(2)
	fsubrl	8(%edi,%esi)
	movl	-1312(%ebp), %esi
	fld	%st(1)
	fldl	8(%eax)
	fld	%st(2)
	fxch	%st(6)
	addl	%esi, %edi
	fstpl	-1328(%ebp)
	fxch	%st(2)
	fmul	%st(4), %st
	fxch	%st(3)
	fmul	%st(2), %st
	fxch	%st(5)
	fmulp	%st, %st(2)
	fmulp	%st, %st(3)
	fxch	%st(3)
	faddp	%st, %st(1)
	fldl	_I
	fldl	_I+8
	fxch	%st(3)
	fsubp	%st, %st(4)
	fld	%st(0)
	fldz
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(3)
	fmulp	%st, %st(4)
	fsubp	%st, %st(2)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(3)
	fadd	%st(1), %st
	fxch	%st(1)
	fstl	-864(%ebp)
	fstpl	-848(%ebp)
	fxch	%st(2)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-880(%ebp)
	fstl	-832(%ebp)
	fxch	%st(1)
	fstl	-856(%ebp)
	fstl	-840(%ebp)
	fstl	-872(%ebp)
	fstl	-824(%ebp)
	fxch	%st(1)
	fsubrl	(%edi)
	fxch	%st(1)
	fsubrl	8(%edi)
	fxch	%st(1)
	fstpl	(%edi)
	fstpl	8(%edi)
L606:
	movl	-1284(%ebp), %esi
	movl	32(%ebp), %eax
	movl	(%eax,%esi,4), %edi
	movl	%edi, -1400(%ebp)
L533:
	movl	-1400(%ebp), %edx
	testl	%edx, %edx
	js	L647
L655:
	incl	-1284(%ebp)
	movl	40(%ebp), %ebx
L654:
	cmpl	%ebx, -1284(%ebp)
	jge	L583
	fldl	8(%ebp)
	fldl	16(%ebp)
	jmp	L581
L425:
	movl	-1404(%ebp), %eax
	movl	-1400(%ebp), %edx
	movl	56(%ebp), %edi
	sall	$4, %eax
	addl	$16, %eax
	imull	%edx, %eax
	addl	%edi, %eax
	xorl	%edi, %edi
	movl	%eax, -1388(%ebp)
	jmp	L638
L648:
	movl	-1300(%ebp), %ecx
	movl	%edi, %eax
	sall	$4, %eax
	fldl	(%eax,%ebx)
	fldl	(%eax,%ecx)
	addl	%eax, %ecx
	fld	%st(1)
	fldl	8(%ecx)
	fxch	%st(1)
	movl	%ecx, -1416(%ebp)
	fmul	%st(2), %st
	movl	48(%ebp), %ecx
	fldl	8(%eax,%ebx)
	fxch	%st(4)
	fmul	%st(2), %st
	leal	(%eax,%ecx), %edx
	fld	%st(4)
	fmulp	%st, %st(4)
	fxch	%st(4)
	fmulp	%st, %st(2)
	fldl	_I+8
	fxch	%st(4)
	faddp	%st, %st(3)
	fldl	_I
	fxch	%st(1)
	fsubp	%st, %st(2)
	fldz
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(4)
	fmulp	%st, %st(5)
	fsubrp	%st, %st(3)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstl	-928(%ebp)
	fstpl	-912(%ebp)
	faddp	%st, %st(2)
	fstl	-944(%ebp)
	fstl	-896(%ebp)
	fxch	%st(1)
	fstl	-920(%ebp)
	fstl	-904(%ebp)
	fstl	-936(%ebp)
	fstl	-888(%ebp)
	fxch	%st(1)
	faddl	(%edx)
	fxch	%st(1)
	faddl	8(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fstpl	8(%edx)
	movl	-1416(%ebp), %edx
	fldl	(%eax,%esi)
	fldl	8(%eax,%esi)
	fldl	(%edx)
	fld	%st(2)
	fld	%st(2)
	fldl	8(%edx)
	fxch	%st(4)
	fmul	%st(3), %st
	fxch	%st(5)
	movl	-1388(%ebp), %edx
	fmul	%st(4), %st
	fxch	%st(1)
	addl	%edx, %eax
	fmulp	%st, %st(4)
	fxch	%st(1)
	fmulp	%st, %st(2)
	faddp	%st, %st(3)
	fldl	_I
	fldl	_I+8
	fxch	%st(2)
	fsubp	%st, %st(3)
	fld	%st(0)
	fldz
	fmul	%st(3), %st
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(5)
	fmulp	%st, %st(3)
	fsubrp	%st, %st(4)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fstl	-992(%ebp)
	fstpl	-976(%ebp)
	faddp	%st, %st(1)
	fstl	-984(%ebp)
	fstl	-968(%ebp)
	fstl	-1000(%ebp)
	fxch	%st(1)
	movl	-1300(%ebp), %edx
	fstl	-1008(%ebp)
	fstl	-960(%ebp)
	fxch	%st(1)
	fstl	-952(%ebp)
	fxch	%st(1)
	faddl	(%eax)
	fxch	%st(1)
	faddl	8(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fstpl	8(%eax)
	leal	1(%edi), %eax
	sall	$4, %eax
	fldl	(%eax,%edx)
	addl	%eax, %edx
	addl	$2, %edi
	fldl	(%eax,%ebx)
	fldl	8(%eax,%ebx)
	fldl	8(%edx)
	fld	%st(2)
	fld	%st(2)
	movl	%edx, -1416(%ebp)
	fmul	%st(5), %st
	fxch	%st(4)
	leal	(%eax,%ecx), %edx
	fmul	%st(2), %st
	fxch	%st(3)
	movl	-1388(%ebp), %ecx
	fmulp	%st, %st(2)
	fmulp	%st, %st(4)
	fxch	%st(1)
	faddp	%st, %st(2)
	fldl	_I
	fldl	_I+8
	fxch	%st(4)
	fsubp	%st, %st(2)
	fld	%st(0)
	fldz
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(4)
	fmulp	%st, %st(5)
	fsubrp	%st, %st(3)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstl	-1056(%ebp)
	fstpl	-1040(%ebp)
	faddp	%st, %st(2)
	fstl	-1072(%ebp)
	fstl	-1024(%ebp)
	fxch	%st(1)
	fstl	-1048(%ebp)
	fstl	-1032(%ebp)
	fstl	-1064(%ebp)
	fstl	-1016(%ebp)
	fxch	%st(1)
	faddl	(%edx)
	fxch	%st(1)
	faddl	8(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fstpl	8(%edx)
	movl	-1416(%ebp), %edx
	fldl	(%eax,%esi)
	fldl	8(%eax,%esi)
	addl	%ecx, %eax
	fld	%st(1)
	fldl	(%edx)
	fld	%st(2)
	fldl	8(%edx)
	fxch	%st(4)
	fmul	%st(2), %st
	fxch	%st(5)
	fmul	%st(4), %st
	fxch	%st(1)
	fmulp	%st, %st(4)
	fxch	%st(2)
	fmulp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(3)
	fldl	_I
	fldl	_I+8
	fxch	%st(2)
	fsubp	%st, %st(3)
	fld	%st(0)
	fldz
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(5)
	fmulp	%st, %st(3)
	fsubp	%st, %st(4)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fstl	-1120(%ebp)
	fstpl	-1104(%ebp)
	faddp	%st, %st(1)
	fxch	%st(1)
	fstl	-1136(%ebp)
	fstl	-1088(%ebp)
	fxch	%st(1)
	fstl	-1112(%ebp)
	fstl	-1096(%ebp)
	fstl	-1128(%ebp)
	fstl	-1080(%ebp)
	fxch	%st(1)
	fsubrl	(%eax)
	fxch	%st(1)
	fsubrl	8(%eax)
	fxch	%st(1)
	fstpl	(%eax)
	fstpl	8(%eax)
L638:
	cmpl	-1404(%ebp), %edi
	jl	L648
	testb	$1, -1404(%ebp)
	jne	L533
	movl	-1404(%ebp), %edi
	movl	-1300(%ebp), %ecx
	sall	$4, %edi
	fldl	(%edi,%ebx)
	leal	(%edi,%ecx), %eax
	fldl	(%edi,%ecx)
	fld	%st(1)
	fldl	8(%edi,%ebx)
	fxch	%st(1)
	fmul	%st(2), %st
	movl	48(%ebp), %ebx
	fldl	8(%eax)
	fld	%st(2)
	fmulp	%st, %st(4)
	leal	(%edi,%ebx), %edx
	fmul	%st, %st(4)
	fmulp	%st, %st(2)
	fldl	_I+8
	fxch	%st(4)
	faddp	%st, %st(3)
	fldl	_I
	fxch	%st(1)
	fsubp	%st, %st(2)
	fldz
	fld	%st(1)
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(4), %st
	fxch	%st(4)
	fmulp	%st, %st(5)
	fsubrp	%st, %st(3)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(1)
	fadd	%st(2), %st
	fxch	%st(2)
	fstl	-1184(%ebp)
	fstpl	-1168(%ebp)
	faddp	%st, %st(2)
	fstl	-1200(%ebp)
	fstl	-1152(%ebp)
	fxch	%st(1)
	fstl	-1176(%ebp)
	fstl	-1160(%ebp)
	fstl	-1192(%ebp)
	fstl	-1144(%ebp)
	fxch	%st(1)
	faddl	(%edx)
	fxch	%st(1)
	faddl	8(%edx)
	fxch	%st(1)
	fstpl	(%edx)
	fstpl	8(%edx)
	fldl	(%edi,%esi)
	fldl	8(%edi,%esi)
	movl	-1388(%ebp), %esi
	fld	%st(1)
	fldl	(%eax)
	fld	%st(2)
	addl	%esi, %edi
	fldl	8(%eax)
	fxch	%st(4)
	fmul	%st(2), %st
	fxch	%st(5)
	fmul	%st(4), %st
	fxch	%st(1)
	fmulp	%st, %st(4)
	fxch	%st(2)
	fmulp	%st, %st(1)
	fxch	%st(1)
	faddp	%st, %st(3)
	fldl	_I
	fldl	_I+8
	fxch	%st(2)
	fsubp	%st, %st(3)
	fld	%st(0)
	fldz
	fxch	%st(1)
	fmul	%st(5), %st
	fxch	%st(1)
	fmul	%st(3), %st
	fxch	%st(5)
	fmulp	%st, %st(3)
	fsubp	%st, %st(4)
	fldz
	fmulp	%st, %st(1)
	fxch	%st(2)
	fadd	%st(3), %st
	fxch	%st(3)
	fstl	-1248(%ebp)
	fstpl	-1232(%ebp)
	faddp	%st, %st(1)
	fstl	-1240(%ebp)
	fstl	-1224(%ebp)
	fstl	-1256(%ebp)
	fxch	%st(1)
	fstl	-1264(%ebp)
	movl	-1400(%ebp), %edx
	fstl	-1216(%ebp)
	fxch	%st(1)
	fstl	-1208(%ebp)
	fxch	%st(1)
	testl	%edx, %edx
	faddl	(%edi)
	fxch	%st(1)
	faddl	8(%edi)
	fxch	%st(1)
	fstpl	(%edi)
	fstpl	8(%edi)
	jns	L655
L647:
	addl	$3, -1284(%ebp)
	movl	40(%ebp), %ebx
	incl	-1284(%ebp)
	jmp	L654
L650:
	fstp	%st(0)
	fstp	%st(0)
	fstp	%st(0)
	jmp	L362
L652:
	fstp	%st(0)
	fstp	%st(1)
	fstp	%st(1)
	jmp	L443
L438:
	leal	-1272(%ebp), %eax
	jmp	L440
L649:
	fstp	%st(0)
	fstp	%st(0)
L583:
	addl	$1420, %esp
	popl	%ebx
	popl	%esi
	popl	%edi
	popl	%ebp
	ret
	.align 2
	.p2align 4,,15
	.def	__GLOBAL__I__Z11shift_allociii;	.scl	3;	.type	32;	.endef
__GLOBAL__I__Z11shift_allociii:
	pushl	%ebp
	fldz
	movl	%esp, %ebp
	popl	%ebp
	fstpl	_I
	fld1
	fstpl	_I+8
	ret
	.def	_atan2;	.scl	2;	.type	32;	.endef
	.def	___fpclassify;	.scl	2;	.type	32;	.endef
	.def	_log;	.scl	2;	.type	32;	.endef
	.def	_mexAtExit;	.scl	2;	.type	32;	.endef
	.def	_exp2;	.scl	2;	.type	32;	.endef
	.def	_mexMakeMemoryPersistent;	.scl	2;	.type	32;	.endef
	.def	_mxMalloc;	.scl	2;	.type	32;	.endef
	.def	_mxFree;	.scl	2;	.type	32;	.endef
