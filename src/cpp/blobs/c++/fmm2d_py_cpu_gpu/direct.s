	.file	"direct.cpp"
	.intel_syntax
	.section	.ctors,"w"
	.align 4
	.long	__GLOBAL__I__Z14directInteractiPKdS0_S0_S0_iS0_S0_PdS1_S1_S1_PK5panelii8SMOOTHERddbS1_i
.lcomm _I,16
.lcomm _TIME_after,16
.lcomm _TIME_before,16
	.section .rdata,"dr"
	.align 4
LC1:
	.long	1056964608
	.align 8
LC2:
	.long	-1670671864
	.long	1072962135
	.align 4
LC3:
	.long	1073741824
	.align 8
LC4:
	.long	-1198367923
	.long	1056062948
	.align 8
LC5:
	.long	1809602294
	.long	1053624430
	.align 8
LC6:
	.long	311951373
	.long	1060283081
	.align 8
LC7:
	.long	-1738515375
	.long	1058950427
	.align 8
LC8:
	.long	272672780
	.long	1064613583
	.align 8
LC9:
	.long	-1945159566
	.long	1063001343
	.align 8
LC10:
	.long	696367074
	.long	1065795916
	.align 8
LC11:
	.long	1293070641
	.long	1066051167
	.align 8
LC12:
	.long	1711857517
	.long	1068693871
	.align 8
LC13:
	.long	-741137424
	.long	1067599537
	.align 8
LC14:
	.long	-321923180
	.long	1072418090
	.align 8
LC15:
	.long	1653800015
	.long	1069699608
	.align 8
LC16:
	.long	-1164543770
	.long	1071747788
	.align 8
LC17:
	.long	2138387848
	.long	1072591351
	.align 8
LC18:
	.long	-1928708653
	.long	1067106279
	.align 8
LC19:
	.long	-705000953
	.long	1072054923
	.align 8
LC20:
	.long	74037635
	.long	1073696940
	.align 8
LC21:
	.long	864351748
	.long	1075789835
	.align 8
LC22:
	.long	-978864670
	.long	1076262673
	.align 8
LC23:
	.long	2102521379
	.long	1077434572
	.align 8
LC24:
	.long	-1419447362
	.long	1077131769
	.align 8
LC25:
	.long	775012524
	.long	1077871833
	.align 8
LC26:
	.long	-231816214
	.long	1076692691
	.align 8
LC27:
	.long	-2117920099
	.long	1077077670
	.align 8
LC28:
	.long	1435360139
	.long	1075073300
	.align 8
LC29:
	.long	-984962297
	.long	1075286286
	.align 8
LC30:
	.long	1307525938
	.long	1072435749
	.align 8
LC31:
	.long	-11994805
	.long	1072580455
	.align 8
LC32:
	.long	275966079
	.long	1068643801
	.align 8
LC33:
	.long	-1776736009
	.long	1068679580
	.align 8
LC34:
	.long	-1261952850
	.long	1063352432
	.align 8
LC35:
	.long	-1261952848
	.long	1063352432
	.align 8
LC36:
	.long	-59787751
	.long	1071806604
	.align 4
LC38:
	.long	-1090519040
	.text
	.align 2
	.p2align 4,,15
.globl __Z5rlog2iPKdS0_S0_PdS1_8SMOOTHERddb
	.def	__Z5rlog2iPKdS0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z5rlog2iPKdS0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 140
	mov	eax, DWORD PTR [ebp+32]
	movzx	edx, BYTE PTR [ebp+52]
	fstp	QWORD PTR [ebp-40]
	mov	edi, DWORD PTR [ebp+24]
	cmp	eax, 1
	je	L12
	jle	L105
	cmp	eax, 2
	je	L37
	cmp	eax, 3
	je	L106
	.p2align 4,,15
L1:
	add	esp, 140
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L37:
	fld	QWORD PTR [ebp+36]
	test	dl, dl
	fmul	st, st(0)
	fst	QWORD PTR [ebp+36]
	fdivr	QWORD PTR LC2
	fstp	QWORD PTR [ebp-96]
	jne	L107
L38:
	xor	esi, esi
	.p2align 4,,15
L101:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L1
	fld	QWORD PTR [ebp-96]
	mov	ebx, DWORD PTR [ebp+20]
	fld	QWORD PTR [ebx+esi*8]
	fmul	DWORD PTR LC1
	fxch	st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fadd	QWORD PTR LC36
	mov	eax, DWORD PTR [ebp+8]
	lea	edx, [esi+1]
	fld	QWORD PTR [ebp-120]
	mov	DWORD PTR [ebp-80], edx
	mov	ebx, edx
	cmp	edx, eax
	fmulp	st(1), st
	fmul	QWORD PTR [ebp-40]
	fadd	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	jge	L76
	.p2align 4,,15
L120:
	mov	eax, DWORD PTR [ebp+16]
	mov	ecx, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+ebx*8]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [ecx+ebx*8]
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+44]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L108
	fldz
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jp	L94
	je	L52
L94:
	fst	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	fld	QWORD PTR [ebp-96]
	fxch	st(2)
	fstp	QWORD PTR [ebp-64]
	fld	DWORD PTR LC3
	fxch	st(2)
	fmul	st, st(1)
	fmulp	st(1), st
	fxch	st(1)
	fucomip	st, st(1)
	jb	L54
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC4
	fxch	st(1)
	fmul	QWORD PTR LC5
	fxch	st(1)
	fadd	QWORD PTR LC6
	fxch	st(1)
	fadd	QWORD PTR LC7
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC8
	fxch	st(1)
	fadd	QWORD PTR LC9
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC10
	fxch	st(1)
	fadd	QWORD PTR LC11
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC12
	fxch	st(1)
	fadd	QWORD PTR LC13
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC14
	fxch	st(1)
	fsub	QWORD PTR LC15
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC16
	fxch	st(1)
	fsub	QWORD PTR LC17
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	fsubrp	st(1), st
L56:
	fmul	DWORD PTR LC1
	fadd	QWORD PTR [ebp-64]
L103:
	fmul	QWORD PTR [ebp-40]
L51:
	mov	edx, DWORD PTR [ebp+20]
	fld	QWORD PTR [edx+ebx*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [edx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L120
L76:
	mov	esi, DWORD PTR [ebp-80]
	jmp	L101
L107:
	fld	QWORD PTR [ebp+44]
	fstp	QWORD PTR [esp]
	call	_log
	fstp	QWORD PTR [ebp-48]
	fld	QWORD PTR [ebp+44]
	fstp	QWORD PTR [esp]
	call	_log
	fld	QWORD PTR [ebp-96]
	fmul	QWORD PTR [ebp+44]
	fxch	st(1)
	fstp	QWORD PTR [ebp-56]
	fld	DWORD PTR LC3
	fxch	st(1)
	fmul	QWORD PTR [ebp+44]
	fxch	st(1)
	fucomip	st, st(1)
	jae	L110
	fld	QWORD PTR [ebp-40]
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC18
	fld	st(1)
	fmul	QWORD PTR LC19
	fxch	st(1)
	fadd	QWORD PTR LC20
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fadd	QWORD PTR LC21
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC22
	fxch	st(2)
	fadd	QWORD PTR LC23
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC24
	fxch	st(2)
	fadd	QWORD PTR LC25
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC26
	fxch	st(2)
	fadd	QWORD PTR LC27
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC28
	fxch	st(2)
	fadd	QWORD PTR LC29
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC30
	fxch	st(2)
	fadd	QWORD PTR LC31
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC32
	fxch	st(2)
	fadd	QWORD PTR LC33
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC34
	fxch	st(1)
	fadd	QWORD PTR LC35
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_exp
	fld	QWORD PTR [ebp-120]
	fld	QWORD PTR [ebp-136]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
L41:
	fmul	DWORD PTR LC1
	fadd	QWORD PTR [ebp-56]
	fdivr	QWORD PTR [ebp-48]
	fstp	QWORD PTR [ebp-40]
	jmp	L38
	.p2align 4,,7
L52:
	fstp	st(0)
	fld	QWORD PTR [ebp-96]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	QWORD PTR LC36
	fmul	DWORD PTR LC38
	jmp	L103
	.p2align 4,,7
L108:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L51
	.p2align 4,,7
L54:
	fld1
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC18
	fld	st(1)
	fmul	QWORD PTR LC19
	fxch	st(1)
	fadd	QWORD PTR LC20
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fadd	QWORD PTR LC21
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC22
	fxch	st(2)
	fadd	QWORD PTR LC23
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC24
	fxch	st(2)
	fadd	QWORD PTR LC25
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC26
	fxch	st(2)
	fadd	QWORD PTR LC27
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC28
	fxch	st(2)
	fadd	QWORD PTR LC29
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC30
	fxch	st(2)
	fadd	QWORD PTR LC31
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC32
	fxch	st(2)
	fadd	QWORD PTR LC33
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC34
	fxch	st(1)
	fadd	QWORD PTR LC35
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_exp
	fld	QWORD PTR [ebp-120]
	fld	QWORD PTR [ebp-136]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
	jmp	L56
L105:
	test	eax, eax
	jne	L1
	fld	QWORD PTR [ebp+36]
	test	dl, dl
	fmul	st, st(0)
	fstp	QWORD PTR [ebp-88]
	jne	L111
L25:
	xor	esi, esi
	.p2align 4,,15
L99:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L1
	fld	QWORD PTR [ebp-88]
	mov	ebx, DWORD PTR [ebp+20]
	fld	QWORD PTR [ebx+esi*8]
	fmul	DWORD PTR LC1
	fxch	st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	mov	edx, DWORD PTR [ebp+8]
	lea	ecx, [esi+1]
	mov	DWORD PTR [ebp-68], ecx
	mov	ebx, ecx
	fmulp	st(1), st
	cmp	ecx, edx
	fmul	QWORD PTR [ebp-40]
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	jge	L73
	.p2align 4,,15
L121:
	mov	edx, DWORD PTR [ebp+16]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [edx+ebx*8]
	fsubr	QWORD PTR [edx+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [eax+ebx*8]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+44]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L32
	fstp	QWORD PTR [esp]
	call	_log
L34:
	mov	ecx, DWORD PTR [ebp+20]
	fld	QWORD PTR [ecx+ebx*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [ecx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L121
L73:
	mov	esi, DWORD PTR [ebp-68]
	jmp	L99
	.p2align 4,,7
L32:
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-88]
	fstp	QWORD PTR [esp]
	call	_log
	fmul	DWORD PTR LC1
	fmul	QWORD PTR [ebp-40]
	jmp	L34
L12:
	fld	QWORD PTR [ebp+36]
	xor	esi, esi
	fstp	QWORD PTR [esp]
	call	_log
	fsub	DWORD PTR LC1
	fstp	QWORD PTR [ebp-32]
	.p2align 4,,15
L97:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L1
	fld	QWORD PTR [ebp-32]
	mov	ebx, DWORD PTR [ebp+20]
	lea	eax, [esi+1]
	mov	DWORD PTR [ebp-76], eax
	mov	ecx, DWORD PTR [ebp+8]
	fmul	QWORD PTR [ebx+esi*8]
	cmp	eax, ecx
	mov	ebx, eax
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	jge	L69
	.p2align 4,,15
L122:
	mov	ecx, DWORD PTR [ebp+16]
	mov	edx, DWORD PTR [ebp+12]
	fld	QWORD PTR [ecx+ebx*8]
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [edx+ebx*8]
	fsubr	QWORD PTR [edx+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+36]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L113
	fdiv	QWORD PTR [ebp+36]
	fld	st(0)
	fmul	DWORD PTR LC1
	fmulp	st(1), st
	fadd	QWORD PTR [ebp-32]
L21:
	mov	eax, DWORD PTR [ebp+20]
	fld	QWORD PTR [eax+ebx*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [eax+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L122
L69:
	mov	esi, DWORD PTR [ebp-76]
	jmp	L97
	.p2align 4,,7
L113:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L21
L106:
	xor	esi, esi
	.p2align 4,,15
L95:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L1
	mov	eax, DWORD PTR [ebp+8]
	lea	edx, [esi+1]
	mov	ebx, edx
	mov	DWORD PTR [ebp-72], edx
	cmp	edx, eax
	jge	L65
	.p2align 4,,15
L123:
	mov	edx, DWORD PTR [ebp+16]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [edx+ebx*8]
	fsubr	QWORD PTR [edx+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [eax+ebx*8]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fstp	QWORD PTR [esp]
	call	_log
	mov	ecx, DWORD PTR [ebp+20]
	fld	QWORD PTR [ecx+ebx*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [ecx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fsubr	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L123
L65:
	mov	esi, DWORD PTR [ebp-72]
	jmp	L95
L110:
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC4
	fxch	st(1)
	fmul	QWORD PTR LC5
	fxch	st(1)
	fadd	QWORD PTR LC6
	fxch	st(1)
	fadd	QWORD PTR LC7
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC8
	fxch	st(1)
	fadd	QWORD PTR LC9
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC10
	fxch	st(1)
	fadd	QWORD PTR LC11
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC12
	fxch	st(1)
	fadd	QWORD PTR LC13
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC14
	fxch	st(1)
	fsub	QWORD PTR LC15
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC16
	fxch	st(1)
	fsub	QWORD PTR LC17
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	fsubrp	st(1), st
	jmp	L41
L111:
	fld	QWORD PTR [ebp+44]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	st, st(0)
	fstp	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp+44]
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-88]
	fstp	QWORD PTR [esp]
	call	_log
	fdivr	QWORD PTR [ebp-40]
	fstp	QWORD PTR [ebp-40]
	jmp	L25
	.section .rdata,"dr"
	.align 4
LC43:
	.long	1056964608
	.align 8
LC44:
	.long	-1670671864
	.long	1072962135
	.align 4
LC45:
	.long	1073741824
	.align 8
LC46:
	.long	-1198367923
	.long	1056062948
	.align 8
LC47:
	.long	1809602294
	.long	1053624430
	.align 8
LC48:
	.long	311951373
	.long	1060283081
	.align 8
LC49:
	.long	-1738515375
	.long	1058950427
	.align 8
LC50:
	.long	272672780
	.long	1064613583
	.align 8
LC51:
	.long	-1945159566
	.long	1063001343
	.align 8
LC52:
	.long	696367074
	.long	1065795916
	.align 8
LC53:
	.long	1293070641
	.long	1066051167
	.align 8
LC54:
	.long	1711857517
	.long	1068693871
	.align 8
LC55:
	.long	-741137424
	.long	1067599537
	.align 8
LC56:
	.long	-321923180
	.long	1072418090
	.align 8
LC57:
	.long	1653800015
	.long	1069699608
	.align 8
LC58:
	.long	-1164543770
	.long	1071747788
	.align 8
LC59:
	.long	2138387848
	.long	1072591351
	.align 8
LC60:
	.long	-1928708653
	.long	1067106279
	.align 8
LC61:
	.long	-705000953
	.long	1072054923
	.align 8
LC62:
	.long	74037635
	.long	1073696940
	.align 8
LC63:
	.long	864351748
	.long	1075789835
	.align 8
LC64:
	.long	-978864670
	.long	1076262673
	.align 8
LC65:
	.long	2102521379
	.long	1077434572
	.align 8
LC66:
	.long	-1419447362
	.long	1077131769
	.align 8
LC67:
	.long	775012524
	.long	1077871833
	.align 8
LC68:
	.long	-231816214
	.long	1076692691
	.align 8
LC69:
	.long	-2117920099
	.long	1077077670
	.align 8
LC70:
	.long	1435360139
	.long	1075073300
	.align 8
LC71:
	.long	-984962297
	.long	1075286286
	.align 8
LC72:
	.long	1307525938
	.long	1072435749
	.align 8
LC73:
	.long	-11994805
	.long	1072580455
	.align 8
LC74:
	.long	275966079
	.long	1068643801
	.align 8
LC75:
	.long	-1776736009
	.long	1068679580
	.align 8
LC76:
	.long	-1261952850
	.long	1063352432
	.align 8
LC77:
	.long	-1261952848
	.long	1063352432
	.align 8
LC78:
	.long	-59787751
	.long	1071806604
	.align 4
LC80:
	.long	-1090519040
	.align 8
LC83:
	.long	0
	.long	1071644672
	.text
	.align 2
	.p2align 4,,15
.globl __Z5zlog2iPKdS0_S0_S0_PdS1_8SMOOTHERddb
	.def	__Z5zlog2iPKdS0_S0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z5zlog2iPKdS0_S0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 140
	mov	eax, DWORD PTR [ebp+36]
	movzx	edx, BYTE PTR [ebp+56]
	fstp	QWORD PTR [ebp-40]
	mov	edi, DWORD PTR [ebp+32]
	cmp	eax, 1
	je	L135
	jle	L225
	cmp	eax, 2
	je	L160
	cmp	eax, 3
	je	L226
	.p2align 4,,15
L124:
	add	esp, 140
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L160:
	fld	QWORD PTR [ebp+40]
	test	dl, dl
	fmul	st, st(0)
	fst	QWORD PTR [ebp+40]
	fdivr	QWORD PTR LC44
	fstp	QWORD PTR [ebp-96]
	jne	L227
L161:
	xor	esi, esi
	.p2align 4,,15
L222:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L124
	fld	QWORD PTR [ebp-96]
	mov	edx, DWORD PTR [ebp+20]
	fld	QWORD PTR [edx+esi*8]
	fmul	QWORD PTR LC83
	fxch	st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fadd	QWORD PTR LC78
	mov	eax, DWORD PTR [ebp+28]
	mov	ebx, DWORD PTR [ebp+24]
	fld	QWORD PTR [ebp-120]
	fmulp	st(1), st
	fld	QWORD PTR [ebp-96]
	fxch	st(1)
	fmul	QWORD PTR [ebp-40]
	fadd	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [eax+esi*8]
	fld	QWORD PTR [ebx+esi*8]
	fmul	QWORD PTR LC83
	fxch	st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fadd	QWORD PTR LC78
	mov	edx, DWORD PTR [ebp+8]
	lea	ecx, [esi+1]
	fld	QWORD PTR [ebp-120]
	mov	DWORD PTR [ebp-76], ecx
	mov	ebx, ecx
	cmp	ecx, edx
	fmulp	st(1), st
	fmul	QWORD PTR [ebp-40]
	fadd	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	jge	L199
	.p2align 4,,15
L240:
	mov	eax, DWORD PTR [ebp+16]
	mov	ecx, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+ebx*8]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [ecx+ebx*8]
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+48]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L228
	fldz
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jp	L215
	je	L175
L215:
	fst	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	fld	QWORD PTR [ebp-96]
	fxch	st(2)
	fstp	QWORD PTR [ebp-64]
	fld	DWORD PTR LC45
	fxch	st(2)
	fmul	st, st(1)
	fmulp	st(1), st
	fxch	st(1)
	fucomip	st, st(1)
	jb	L177
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC46
	fxch	st(1)
	fmul	QWORD PTR LC47
	fxch	st(1)
	fadd	QWORD PTR LC48
	fxch	st(1)
	fadd	QWORD PTR LC49
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC50
	fxch	st(1)
	fadd	QWORD PTR LC51
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC52
	fxch	st(1)
	fadd	QWORD PTR LC53
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC54
	fxch	st(1)
	fadd	QWORD PTR LC55
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC56
	fxch	st(1)
	fsub	QWORD PTR LC57
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC58
	fxch	st(1)
	fsub	QWORD PTR LC59
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	fsubrp	st(1), st
L179:
	fmul	DWORD PTR LC43
	fadd	QWORD PTR [ebp-64]
L224:
	fmul	QWORD PTR [ebp-40]
L174:
	mov	eax, DWORD PTR [ebp+20]
	mov	ecx, DWORD PTR [ebp+28]
	mov	edx, DWORD PTR [ebp+24]
	fld	QWORD PTR [eax+ebx*8]
	mov	eax, DWORD PTR [ebp+20]
	fmul	st, st(1)
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [ecx+esi*8]
	mov	ecx, DWORD PTR [ebp+28]
	fld	QWORD PTR [edx+ebx*8]
	mov	edx, DWORD PTR [ebp+24]
	fmul	st, st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fld	QWORD PTR [eax+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [edx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [ecx+ebx*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L240
L199:
	mov	esi, DWORD PTR [ebp-76]
	jmp	L222
L227:
	fld	QWORD PTR [ebp+48]
	fstp	QWORD PTR [esp]
	call	_log
	fstp	QWORD PTR [ebp-48]
	fld	QWORD PTR [ebp+48]
	fstp	QWORD PTR [esp]
	call	_log
	fld	QWORD PTR [ebp-96]
	fmul	QWORD PTR [ebp+48]
	fxch	st(1)
	fstp	QWORD PTR [ebp-56]
	fld	DWORD PTR LC45
	fxch	st(1)
	fmul	QWORD PTR [ebp+48]
	fxch	st(1)
	fucomip	st, st(1)
	jae	L230
	fld	QWORD PTR [ebp-40]
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC60
	fld	st(1)
	fmul	QWORD PTR LC61
	fxch	st(1)
	fadd	QWORD PTR LC62
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fadd	QWORD PTR LC63
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC64
	fxch	st(2)
	fadd	QWORD PTR LC65
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC66
	fxch	st(2)
	fadd	QWORD PTR LC67
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC68
	fxch	st(2)
	fadd	QWORD PTR LC69
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC70
	fxch	st(2)
	fadd	QWORD PTR LC71
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC72
	fxch	st(2)
	fadd	QWORD PTR LC73
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC74
	fxch	st(2)
	fadd	QWORD PTR LC75
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC76
	fxch	st(1)
	fadd	QWORD PTR LC77
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_exp
	fld	QWORD PTR [ebp-120]
	fld	QWORD PTR [ebp-136]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
L164:
	fmul	DWORD PTR LC43
	fadd	QWORD PTR [ebp-56]
	fdivr	QWORD PTR [ebp-48]
	fstp	QWORD PTR [ebp-40]
	jmp	L161
	.p2align 4,,7
L175:
	fstp	st(0)
	fld	QWORD PTR [ebp-96]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	QWORD PTR LC78
	fmul	DWORD PTR LC80
	jmp	L224
	.p2align 4,,7
L228:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L174
	.p2align 4,,7
L177:
	fld1
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC60
	fld	st(1)
	fmul	QWORD PTR LC61
	fxch	st(1)
	fadd	QWORD PTR LC62
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fadd	QWORD PTR LC63
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC64
	fxch	st(2)
	fadd	QWORD PTR LC65
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC66
	fxch	st(2)
	fadd	QWORD PTR LC67
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC68
	fxch	st(2)
	fadd	QWORD PTR LC69
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC70
	fxch	st(2)
	fadd	QWORD PTR LC71
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC72
	fxch	st(2)
	fadd	QWORD PTR LC73
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC74
	fxch	st(2)
	fadd	QWORD PTR LC75
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC76
	fxch	st(1)
	fadd	QWORD PTR LC77
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_exp
	fld	QWORD PTR [ebp-120]
	fld	QWORD PTR [ebp-136]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
	jmp	L179
L225:
	test	eax, eax
	jne	L124
	fld	QWORD PTR [ebp+40]
	test	dl, dl
	fmul	st, st(0)
	fstp	QWORD PTR [ebp-88]
	jne	L231
L148:
	xor	esi, esi
	.p2align 4,,15
L220:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L124
	fld	QWORD PTR [ebp-88]
	mov	edx, DWORD PTR [ebp+20]
	fld	QWORD PTR [edx+esi*8]
	fmul	QWORD PTR LC83
	fxch	st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	mov	eax, DWORD PTR [ebp+28]
	mov	ebx, DWORD PTR [ebp+24]
	fmulp	st(1), st
	fld	QWORD PTR [ebp-88]
	fxch	st(1)
	fmul	QWORD PTR [ebp-40]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [eax+esi*8]
	fld	QWORD PTR [ebx+esi*8]
	fmul	QWORD PTR LC83
	fxch	st(1)
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	mov	edx, DWORD PTR [ebp+8]
	lea	ecx, [esi+1]
	mov	DWORD PTR [ebp-80], ecx
	mov	ebx, ecx
	fmulp	st(1), st
	cmp	ecx, edx
	fmul	QWORD PTR [ebp-40]
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	jge	L196
	.p2align 4,,15
L241:
	mov	eax, DWORD PTR [ebp+16]
	mov	ecx, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+ebx*8]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [ecx+ebx*8]
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+48]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L155
	fstp	QWORD PTR [esp]
	call	_log
L157:
	mov	eax, DWORD PTR [ebp+20]
	mov	ecx, DWORD PTR [ebp+28]
	mov	edx, DWORD PTR [ebp+24]
	fld	QWORD PTR [eax+ebx*8]
	mov	eax, DWORD PTR [ebp+20]
	fmul	st, st(1)
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [ecx+esi*8]
	mov	ecx, DWORD PTR [ebp+28]
	fld	QWORD PTR [edx+ebx*8]
	mov	edx, DWORD PTR [ebp+24]
	fmul	st, st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fld	QWORD PTR [eax+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [edx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [ecx+ebx*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L241
L196:
	mov	esi, DWORD PTR [ebp-80]
	jmp	L220
	.p2align 4,,7
L155:
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-88]
	fstp	QWORD PTR [esp]
	call	_log
	fmul	DWORD PTR LC43
	fmul	QWORD PTR [ebp-40]
	jmp	L157
L135:
	fld	QWORD PTR [ebp+40]
	xor	esi, esi
	fstp	QWORD PTR [esp]
	call	_log
	fsub	DWORD PTR LC43
	fstp	QWORD PTR [ebp-32]
	.p2align 4,,15
L218:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L124
	fld	QWORD PTR [ebp-32]
	mov	edx, DWORD PTR [ebp+20]
	lea	ecx, [esi+1]
	mov	DWORD PTR [ebp-72], ecx
	mov	eax, DWORD PTR [ebp+28]
	mov	ebx, DWORD PTR [ebp+24]
	fmul	QWORD PTR [edx+esi*8]
	mov	edx, DWORD PTR [ebp+8]
	cmp	ecx, edx
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [eax+esi*8]
	fld	QWORD PTR [ebp-32]
	fmul	QWORD PTR [ebx+esi*8]
	mov	ebx, ecx
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	jge	L192
	.p2align 4,,15
L242:
	mov	eax, DWORD PTR [ebp+16]
	mov	ecx, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+ebx*8]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [ecx+ebx*8]
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+40]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L233
	fdiv	QWORD PTR [ebp+40]
	fld	st(0)
	fmul	DWORD PTR LC43
	fmulp	st(1), st
	fadd	QWORD PTR [ebp-32]
L144:
	mov	eax, DWORD PTR [ebp+20]
	mov	ecx, DWORD PTR [ebp+28]
	mov	edx, DWORD PTR [ebp+24]
	fld	QWORD PTR [eax+ebx*8]
	mov	eax, DWORD PTR [ebp+20]
	fmul	st, st(1)
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [ecx+esi*8]
	mov	ecx, DWORD PTR [ebp+28]
	fld	QWORD PTR [edx+ebx*8]
	mov	edx, DWORD PTR [ebp+24]
	fmul	st, st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fld	QWORD PTR [eax+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [edx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [ecx+ebx*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L242
L192:
	mov	esi, DWORD PTR [ebp-72]
	jmp	L218
	.p2align 4,,7
L233:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L144
L226:
	xor	esi, esi
	.p2align 4,,15
L216:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L124
	mov	eax, DWORD PTR [ebp+8]
	lea	edx, [esi+1]
	mov	ebx, edx
	mov	DWORD PTR [ebp-68], edx
	cmp	edx, eax
	jge	L188
	.p2align 4,,15
L243:
	mov	eax, DWORD PTR [ebp+16]
	mov	ecx, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+ebx*8]
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [ecx+ebx*8]
	fsubr	QWORD PTR [ecx+esi*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fstp	QWORD PTR [esp]
	call	_log
	mov	edx, DWORD PTR [ebp+20]
	mov	eax, DWORD PTR [ebp+28]
	mov	ecx, DWORD PTR [ebp+24]
	fld	QWORD PTR [edx+ebx*8]
	mov	edx, DWORD PTR [ebp+20]
	fmul	st, st(1)
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [eax+esi*8]
	mov	eax, DWORD PTR [ebp+28]
	fld	QWORD PTR [ecx+ebx*8]
	mov	ecx, DWORD PTR [ebp+24]
	fmul	st, st(1)
	fsubr	QWORD PTR [edi+esi*8]
	fstp	QWORD PTR [edi+esi*8]
	fld	QWORD PTR [edx+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [ecx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [eax+ebx*8]
	fxch	st(1)
	fsubr	QWORD PTR [edi+ebx*8]
	fxch	st(1)
	fstp	QWORD PTR [eax+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L243
L188:
	mov	esi, DWORD PTR [ebp-68]
	jmp	L216
L230:
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC46
	fxch	st(1)
	fmul	QWORD PTR LC47
	fxch	st(1)
	fadd	QWORD PTR LC48
	fxch	st(1)
	fadd	QWORD PTR LC49
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC50
	fxch	st(1)
	fadd	QWORD PTR LC51
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC52
	fxch	st(1)
	fadd	QWORD PTR LC53
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC54
	fxch	st(1)
	fadd	QWORD PTR LC55
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC56
	fxch	st(1)
	fsub	QWORD PTR LC57
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC58
	fxch	st(1)
	fsub	QWORD PTR LC59
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-120]
	call	_log
	fld	QWORD PTR [ebp-120]
	fsubrp	st(1), st
	jmp	L164
L231:
	fld	QWORD PTR [ebp+48]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	st, st(0)
	fstp	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp+48]
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-88]
	fstp	QWORD PTR [esp]
	call	_log
	fdivr	QWORD PTR [ebp-40]
	fstp	QWORD PTR [ebp-40]
	jmp	L148
	.section .rdata,"dr"
	.align 4
LC86:
	.long	1056964608
	.align 8
LC87:
	.long	-1670671864
	.long	1072962135
	.align 4
LC88:
	.long	1073741824
	.align 8
LC89:
	.long	-1198367923
	.long	1056062948
	.align 8
LC90:
	.long	1809602294
	.long	1053624430
	.align 8
LC91:
	.long	311951373
	.long	1060283081
	.align 8
LC92:
	.long	-1738515375
	.long	1058950427
	.align 8
LC93:
	.long	272672780
	.long	1064613583
	.align 8
LC94:
	.long	-1945159566
	.long	1063001343
	.align 8
LC95:
	.long	696367074
	.long	1065795916
	.align 8
LC96:
	.long	1293070641
	.long	1066051167
	.align 8
LC97:
	.long	1711857517
	.long	1068693871
	.align 8
LC98:
	.long	-741137424
	.long	1067599537
	.align 8
LC99:
	.long	-321923180
	.long	1072418090
	.align 8
LC100:
	.long	1653800015
	.long	1069699608
	.align 8
LC101:
	.long	-1164543770
	.long	1071747788
	.align 8
LC102:
	.long	2138387848
	.long	1072591351
	.align 8
LC103:
	.long	-1928708653
	.long	1067106279
	.align 8
LC104:
	.long	-705000953
	.long	1072054923
	.align 8
LC105:
	.long	74037635
	.long	1073696940
	.align 8
LC106:
	.long	864351748
	.long	1075789835
	.align 8
LC107:
	.long	-978864670
	.long	1076262673
	.align 8
LC108:
	.long	2102521379
	.long	1077434572
	.align 8
LC109:
	.long	-1419447362
	.long	1077131769
	.align 8
LC110:
	.long	775012524
	.long	1077871833
	.align 8
LC111:
	.long	-231816214
	.long	1076692691
	.align 8
LC112:
	.long	-2117920099
	.long	1077077670
	.align 8
LC113:
	.long	1435360139
	.long	1075073300
	.align 8
LC114:
	.long	-984962297
	.long	1075286286
	.align 8
LC115:
	.long	1307525938
	.long	1072435749
	.align 8
LC116:
	.long	-11994805
	.long	1072580455
	.align 8
LC117:
	.long	275966079
	.long	1068643801
	.align 8
LC118:
	.long	-1776736009
	.long	1068679580
	.align 8
LC119:
	.long	-1261952850
	.long	1063352432
	.align 8
LC120:
	.long	-1261952848
	.long	1063352432
	.align 8
LC122:
	.long	-59787751
	.long	1071806604
	.align 4
LC123:
	.long	-1090519040
	.text
	.align 2
	.p2align 4,,15
.globl __Z4rlogiPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb
	.def	__Z4rlogiPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z4rlogiPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 124
	mov	eax, DWORD PTR [ebp+44]
	movzx	edx, BYTE PTR [ebp+64]
	fstp	QWORD PTR [ebp-40]
	mov	edi, DWORD PTR [ebp+36]
	cmp	eax, 1
	je	L257
	jle	L351
	cmp	eax, 2
	je	L282
	cmp	eax, 3
	je	L352
	.p2align 4,,15
L244:
	add	esp, 124
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L282:
	fld	QWORD PTR [ebp+48]
	test	dl, dl
	fmul	st, st(0)
	fst	QWORD PTR [ebp+48]
	fdivr	QWORD PTR LC87
	fstp	QWORD PTR [ebp-80]
	jne	L353
L283:
	xor	esi, esi
	.p2align 4,,15
L348:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L244
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jge	L321
	.p2align 4,,15
L368:
	mov	ecx, DWORD PTR [ebp+16]
	mov	edx, DWORD PTR [ebp+32]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ecx+esi*8]
	mov	ecx, DWORD PTR [ebp+28]
	fsubr	QWORD PTR [edx+ebx*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [eax+esi*8]
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+56]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L354
	fldz
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jp	L341
	je	L297
L341:
	fst	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-104]
	call	_log
	fld	QWORD PTR [ebp-104]
	fld	QWORD PTR [ebp-80]
	fxch	st(2)
	fstp	QWORD PTR [ebp-64]
	fld	DWORD PTR LC88
	fxch	st(2)
	fmul	st, st(1)
	fmulp	st(1), st
	fxch	st(1)
	fucomip	st, st(1)
	jb	L299
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC89
	fxch	st(1)
	fmul	QWORD PTR LC90
	fxch	st(1)
	fadd	QWORD PTR LC91
	fxch	st(1)
	fadd	QWORD PTR LC92
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC93
	fxch	st(1)
	fadd	QWORD PTR LC94
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC95
	fxch	st(1)
	fadd	QWORD PTR LC96
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC97
	fxch	st(1)
	fadd	QWORD PTR LC98
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC99
	fxch	st(1)
	fsub	QWORD PTR LC100
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC101
	fxch	st(1)
	fsub	QWORD PTR LC102
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_log
	fld	QWORD PTR [ebp-104]
	fsubrp	st(1), st
L301:
	fmul	DWORD PTR LC86
	fadd	QWORD PTR [ebp-64]
L350:
	fmul	QWORD PTR [ebp-40]
L296:
	mov	eax, DWORD PTR [ebp+20]
	fmul	QWORD PTR [eax+esi*8]
	fsubr	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jl	L368
L321:
	inc	esi
	jmp	L348
L353:
	fld	QWORD PTR [ebp+56]
	fstp	QWORD PTR [esp]
	call	_log
	fstp	QWORD PTR [ebp-48]
	fld	QWORD PTR [ebp+56]
	fstp	QWORD PTR [esp]
	call	_log
	fld	QWORD PTR [ebp-80]
	fmul	QWORD PTR [ebp+56]
	fxch	st(1)
	fstp	QWORD PTR [ebp-56]
	fld	DWORD PTR LC88
	fxch	st(1)
	fmul	QWORD PTR [ebp+56]
	fxch	st(1)
	fucomip	st, st(1)
	jae	L356
	fld	QWORD PTR [ebp-40]
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC103
	fld	st(1)
	fmul	QWORD PTR LC104
	fxch	st(1)
	fadd	QWORD PTR LC105
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-120]
	fxch	st(1)
	fadd	QWORD PTR LC106
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC107
	fxch	st(2)
	fadd	QWORD PTR LC108
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC109
	fxch	st(2)
	fadd	QWORD PTR LC110
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC111
	fxch	st(2)
	fadd	QWORD PTR LC112
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC113
	fxch	st(2)
	fadd	QWORD PTR LC114
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC115
	fxch	st(2)
	fadd	QWORD PTR LC116
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC117
	fxch	st(2)
	fadd	QWORD PTR LC118
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC119
	fxch	st(1)
	fadd	QWORD PTR LC120
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_exp
	fld	QWORD PTR [ebp-104]
	fld	QWORD PTR [ebp-120]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
L286:
	fmul	DWORD PTR LC86
	fadd	QWORD PTR [ebp-56]
	fdivr	QWORD PTR [ebp-48]
	fstp	QWORD PTR [ebp-40]
	jmp	L283
	.p2align 4,,7
L297:
	fstp	st(0)
	fld	QWORD PTR [ebp-80]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	QWORD PTR LC122
	fmul	DWORD PTR LC123
	jmp	L350
	.p2align 4,,7
L354:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L296
	.p2align 4,,7
L299:
	fld1
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC103
	fld	st(1)
	fmul	QWORD PTR LC104
	fxch	st(1)
	fadd	QWORD PTR LC105
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-120]
	fxch	st(1)
	fadd	QWORD PTR LC106
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC107
	fxch	st(2)
	fadd	QWORD PTR LC108
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC109
	fxch	st(2)
	fadd	QWORD PTR LC110
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC111
	fxch	st(2)
	fadd	QWORD PTR LC112
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC113
	fxch	st(2)
	fadd	QWORD PTR LC114
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC115
	fxch	st(2)
	fadd	QWORD PTR LC116
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC117
	fxch	st(2)
	fadd	QWORD PTR LC118
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC119
	fxch	st(1)
	fadd	QWORD PTR LC120
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_exp
	fld	QWORD PTR [ebp-104]
	fld	QWORD PTR [ebp-120]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
	jmp	L301
L351:
	test	eax, eax
	jne	L244
	fld	QWORD PTR [ebp+48]
	test	dl, dl
	fmul	st, st(0)
	fstp	QWORD PTR [ebp-72]
	jne	L357
L270:
	xor	esi, esi
	.p2align 4,,15
L346:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L244
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jge	L318
	.p2align 4,,15
L369:
	mov	eax, DWORD PTR [ebp+16]
	mov	ecx, DWORD PTR [ebp+32]
	mov	edx, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+esi*8]
	mov	eax, DWORD PTR [ebp+28]
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [edx+esi*8]
	fsubr	QWORD PTR [eax+ebx*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+56]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L277
	fstp	QWORD PTR [esp]
	call	_log
L279:
	mov	edx, DWORD PTR [ebp+20]
	fmul	QWORD PTR [edx+esi*8]
	fsubr	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jl	L369
L318:
	inc	esi
	jmp	L346
	.p2align 4,,7
L277:
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-72]
	fstp	QWORD PTR [esp]
	call	_log
	fmul	DWORD PTR LC86
	fmul	QWORD PTR [ebp-40]
	jmp	L279
L257:
	fld	QWORD PTR [ebp+48]
	xor	esi, esi
	fstp	QWORD PTR [esp]
	call	_log
	fsub	DWORD PTR LC86
	fstp	QWORD PTR [ebp-32]
	.p2align 4,,15
L344:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L244
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jge	L314
	.p2align 4,,15
L370:
	mov	edx, DWORD PTR [ebp+16]
	mov	eax, DWORD PTR [ebp+32]
	mov	ecx, DWORD PTR [ebp+12]
	fld	QWORD PTR [edx+esi*8]
	mov	edx, DWORD PTR [ebp+28]
	fsubr	QWORD PTR [eax+ebx*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [ecx+esi*8]
	fsubr	QWORD PTR [edx+ebx*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+48]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L359
	fdiv	QWORD PTR [ebp+48]
	fld	st(0)
	fmul	DWORD PTR LC86
	fmulp	st(1), st
	fadd	QWORD PTR [ebp-32]
L266:
	fld	QWORD PTR [edi+ebx*8]
	fxch	st(1)
	mov	ecx, DWORD PTR [ebp+20]
	fmul	QWORD PTR [ecx+esi*8]
	fsubp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jl	L370
L314:
	inc	esi
	jmp	L344
	.p2align 4,,7
L359:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L266
L352:
	xor	esi, esi
	.p2align 4,,15
L342:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L244
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jge	L310
	.p2align 4,,15
L371:
	mov	edx, DWORD PTR [ebp+28]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [edx+ebx*8]
	fld	QWORD PTR [eax+esi*8]
	fxch	st(1)
	fucomi	st, st(1)
	jne	L322
	jp	L322
	mov	eax, DWORD PTR [ebp+32]
	mov	ecx, DWORD PTR [ebp+16]
	fld	QWORD PTR [eax+ebx*8]
	fld	QWORD PTR [ecx+esi*8]
	fxch	st(1)
	fucomi	st, st(1)
	jp	L362
	je	L363
L362:
	fxch	st(1)
	fsubp	st(1), st
	fxch	st(1)
	fsubrp	st(2), st
	fstp	QWORD PTR [esp+8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fstp	QWORD PTR [esp]
	call	_log
	mov	eax, DWORD PTR [ebp+20]
	fmul	QWORD PTR [eax+esi*8]
	fsubr	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	.p2align 4,,15
L372:
	cmp	ebx, DWORD PTR [ebp+24]
	jl	L371
L310:
	inc	esi
	jmp	L342
	.p2align 4,,7
L363:
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	inc	ebx
	jmp	L372
	.p2align 4,,7
L322:
	mov	ecx, DWORD PTR [ebp+32]
	mov	edx, DWORD PTR [ebp+16]
	fld	QWORD PTR [ecx+ebx*8]
	fld	QWORD PTR [edx+esi*8]
	fsubp	st(1), st
	fxch	st(1)
	fsubrp	st(2), st
	fstp	QWORD PTR [esp+8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fstp	QWORD PTR [esp]
	call	_log
	mov	eax, DWORD PTR [ebp+20]
	fmul	QWORD PTR [eax+esi*8]
	fsubr	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	jmp	L372
L356:
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC89
	fxch	st(1)
	fmul	QWORD PTR LC90
	fxch	st(1)
	fadd	QWORD PTR LC91
	fxch	st(1)
	fadd	QWORD PTR LC92
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC93
	fxch	st(1)
	fadd	QWORD PTR LC94
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC95
	fxch	st(1)
	fadd	QWORD PTR LC96
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC97
	fxch	st(1)
	fadd	QWORD PTR LC98
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC99
	fxch	st(1)
	fsub	QWORD PTR LC100
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC101
	fxch	st(1)
	fsub	QWORD PTR LC102
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_log
	fld	QWORD PTR [ebp-104]
	fsubrp	st(1), st
	jmp	L286
L357:
	fld	QWORD PTR [ebp+56]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	st, st(0)
	fstp	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp+56]
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-72]
	fstp	QWORD PTR [esp]
	call	_log
	fdivr	QWORD PTR [ebp-40]
	fstp	QWORD PTR [ebp-40]
	jmp	L270
	.section .rdata,"dr"
	.align 4
LC128:
	.long	1056964608
	.align 8
LC129:
	.long	-1670671864
	.long	1072962135
	.align 4
LC130:
	.long	1073741824
	.align 8
LC131:
	.long	-1198367923
	.long	1056062948
	.align 8
LC132:
	.long	1809602294
	.long	1053624430
	.align 8
LC133:
	.long	311951373
	.long	1060283081
	.align 8
LC134:
	.long	-1738515375
	.long	1058950427
	.align 8
LC135:
	.long	272672780
	.long	1064613583
	.align 8
LC136:
	.long	-1945159566
	.long	1063001343
	.align 8
LC137:
	.long	696367074
	.long	1065795916
	.align 8
LC138:
	.long	1293070641
	.long	1066051167
	.align 8
LC139:
	.long	1711857517
	.long	1068693871
	.align 8
LC140:
	.long	-741137424
	.long	1067599537
	.align 8
LC141:
	.long	-321923180
	.long	1072418090
	.align 8
LC142:
	.long	1653800015
	.long	1069699608
	.align 8
LC143:
	.long	-1164543770
	.long	1071747788
	.align 8
LC144:
	.long	2138387848
	.long	1072591351
	.align 8
LC145:
	.long	-1928708653
	.long	1067106279
	.align 8
LC146:
	.long	-705000953
	.long	1072054923
	.align 8
LC147:
	.long	74037635
	.long	1073696940
	.align 8
LC148:
	.long	864351748
	.long	1075789835
	.align 8
LC149:
	.long	-978864670
	.long	1076262673
	.align 8
LC150:
	.long	2102521379
	.long	1077434572
	.align 8
LC151:
	.long	-1419447362
	.long	1077131769
	.align 8
LC152:
	.long	775012524
	.long	1077871833
	.align 8
LC153:
	.long	-231816214
	.long	1076692691
	.align 8
LC154:
	.long	-2117920099
	.long	1077077670
	.align 8
LC155:
	.long	1435360139
	.long	1075073300
	.align 8
LC156:
	.long	-984962297
	.long	1075286286
	.align 8
LC157:
	.long	1307525938
	.long	1072435749
	.align 8
LC158:
	.long	-11994805
	.long	1072580455
	.align 8
LC159:
	.long	275966079
	.long	1068643801
	.align 8
LC160:
	.long	-1776736009
	.long	1068679580
	.align 8
LC161:
	.long	-1261952850
	.long	1063352432
	.align 8
LC162:
	.long	-1261952848
	.long	1063352432
	.align 8
LC164:
	.long	-59787751
	.long	1071806604
	.align 4
LC165:
	.long	-1090519040
	.text
	.align 2
	.p2align 4,,15
.globl __Z4zlogiPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb
	.def	__Z4zlogiPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z4zlogiPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 124
	mov	eax, DWORD PTR [ebp+48]
	movzx	edx, BYTE PTR [ebp+68]
	fstp	QWORD PTR [ebp-40]
	mov	edi, DWORD PTR [ebp+44]
	cmp	eax, 1
	je	L386
	jle	L480
	cmp	eax, 2
	je	L411
	cmp	eax, 3
	je	L481
	.p2align 4,,15
L373:
	add	esp, 124
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L411:
	fld	QWORD PTR [ebp+52]
	test	dl, dl
	fmul	st, st(0)
	fst	QWORD PTR [ebp+52]
	fdivr	QWORD PTR LC129
	fstp	QWORD PTR [ebp-80]
	jne	L482
L412:
	xor	esi, esi
	.p2align 4,,15
L477:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L373
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jge	L450
	.p2align 4,,15
L497:
	mov	ecx, DWORD PTR [ebp+16]
	mov	edx, DWORD PTR [ebp+36]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ecx+esi*8]
	mov	ecx, DWORD PTR [ebp+32]
	fsubr	QWORD PTR [edx+ebx*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [eax+esi*8]
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+60]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L483
	fldz
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jp	L470
	je	L426
L470:
	fst	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-104]
	call	_log
	fld	QWORD PTR [ebp-104]
	fld	QWORD PTR [ebp-80]
	fxch	st(2)
	fstp	QWORD PTR [ebp-64]
	fld	DWORD PTR LC130
	fxch	st(2)
	fmul	st, st(1)
	fmulp	st(1), st
	fxch	st(1)
	fucomip	st, st(1)
	jb	L428
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC131
	fxch	st(1)
	fmul	QWORD PTR LC132
	fxch	st(1)
	fadd	QWORD PTR LC133
	fxch	st(1)
	fadd	QWORD PTR LC134
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC135
	fxch	st(1)
	fadd	QWORD PTR LC136
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC137
	fxch	st(1)
	fadd	QWORD PTR LC138
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC139
	fxch	st(1)
	fadd	QWORD PTR LC140
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC141
	fxch	st(1)
	fsub	QWORD PTR LC142
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC143
	fxch	st(1)
	fsub	QWORD PTR LC144
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_log
	fld	QWORD PTR [ebp-104]
	fsubrp	st(1), st
L430:
	fmul	DWORD PTR LC128
	fadd	QWORD PTR [ebp-64]
L479:
	fmul	QWORD PTR [ebp-40]
L425:
	mov	ecx, DWORD PTR [ebp+20]
	mov	edx, DWORD PTR [ebp+40]
	mov	eax, DWORD PTR [ebp+24]
	fld	QWORD PTR [ecx+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [eax+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [edx+ebx*8]
	fstp	QWORD PTR [edx+ebx*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jl	L497
L450:
	inc	esi
	jmp	L477
L482:
	fld	QWORD PTR [ebp+60]
	fstp	QWORD PTR [esp]
	call	_log
	fstp	QWORD PTR [ebp-48]
	fld	QWORD PTR [ebp+60]
	fstp	QWORD PTR [esp]
	call	_log
	fld	QWORD PTR [ebp-80]
	fmul	QWORD PTR [ebp+60]
	fxch	st(1)
	fstp	QWORD PTR [ebp-56]
	fld	DWORD PTR LC130
	fxch	st(1)
	fmul	QWORD PTR [ebp+60]
	fxch	st(1)
	fucomip	st, st(1)
	jae	L485
	fld	QWORD PTR [ebp-40]
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC145
	fld	st(1)
	fmul	QWORD PTR LC146
	fxch	st(1)
	fadd	QWORD PTR LC147
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-120]
	fxch	st(1)
	fadd	QWORD PTR LC148
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC149
	fxch	st(2)
	fadd	QWORD PTR LC150
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC151
	fxch	st(2)
	fadd	QWORD PTR LC152
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC153
	fxch	st(2)
	fadd	QWORD PTR LC154
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC155
	fxch	st(2)
	fadd	QWORD PTR LC156
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC157
	fxch	st(2)
	fadd	QWORD PTR LC158
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC159
	fxch	st(2)
	fadd	QWORD PTR LC160
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC161
	fxch	st(1)
	fadd	QWORD PTR LC162
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_exp
	fld	QWORD PTR [ebp-104]
	fld	QWORD PTR [ebp-120]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
L415:
	fmul	DWORD PTR LC128
	fadd	QWORD PTR [ebp-56]
	fdivr	QWORD PTR [ebp-48]
	fstp	QWORD PTR [ebp-40]
	jmp	L412
	.p2align 4,,7
L426:
	fstp	st(0)
	fld	QWORD PTR [ebp-80]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	QWORD PTR LC164
	fmul	DWORD PTR LC165
	jmp	L479
	.p2align 4,,7
L483:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L425
	.p2align 4,,7
L428:
	fld1
	fdiv	st, st(1)
	fld	st(0)
	fmul	QWORD PTR LC145
	fld	st(1)
	fmul	QWORD PTR LC146
	fxch	st(1)
	fadd	QWORD PTR LC147
	fxch	st(3)
	fchs
	fxch	st(2)
	fst	QWORD PTR [ebp-120]
	fxch	st(1)
	fadd	QWORD PTR LC148
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(2)
	fstp	QWORD PTR [esp]
	fmul	st(2), st
	fxch	st(1)
	fadd	QWORD PTR LC149
	fxch	st(2)
	fadd	QWORD PTR LC150
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC151
	fxch	st(2)
	fadd	QWORD PTR LC152
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC153
	fxch	st(2)
	fadd	QWORD PTR LC154
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC155
	fxch	st(2)
	fadd	QWORD PTR LC156
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC157
	fxch	st(2)
	fadd	QWORD PTR LC158
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fadd	QWORD PTR LC159
	fxch	st(2)
	fadd	QWORD PTR LC160
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR LC161
	fxch	st(1)
	fadd	QWORD PTR LC162
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_exp
	fld	QWORD PTR [ebp-104]
	fld	QWORD PTR [ebp-120]
	fxch	st(2)
	fmulp	st(1), st
	fmulp	st(1), st
	jmp	L430
L480:
	test	eax, eax
	jne	L373
	fld	QWORD PTR [ebp+52]
	test	dl, dl
	fmul	st, st(0)
	fstp	QWORD PTR [ebp-72]
	jne	L486
L399:
	xor	esi, esi
	.p2align 4,,15
L475:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L373
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jge	L447
	.p2align 4,,15
L498:
	mov	edx, DWORD PTR [ebp+16]
	mov	eax, DWORD PTR [ebp+36]
	mov	ecx, DWORD PTR [ebp+12]
	fld	QWORD PTR [edx+esi*8]
	mov	edx, DWORD PTR [ebp+32]
	fsubr	QWORD PTR [eax+ebx*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [ecx+esi*8]
	fsubr	QWORD PTR [edx+ebx*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+60]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L406
	fstp	QWORD PTR [esp]
	call	_log
L408:
	mov	edx, DWORD PTR [ebp+20]
	mov	eax, DWORD PTR [ebp+40]
	mov	ecx, DWORD PTR [ebp+24]
	fld	QWORD PTR [edx+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [ecx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [eax+ebx*8]
	fstp	QWORD PTR [eax+ebx*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jl	L498
L447:
	inc	esi
	jmp	L475
	.p2align 4,,7
L406:
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-72]
	fstp	QWORD PTR [esp]
	call	_log
	fmul	DWORD PTR LC128
	fmul	QWORD PTR [ebp-40]
	jmp	L408
L386:
	fld	QWORD PTR [ebp+52]
	xor	esi, esi
	fstp	QWORD PTR [esp]
	call	_log
	fsub	DWORD PTR LC128
	fstp	QWORD PTR [ebp-32]
	.p2align 4,,15
L473:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L373
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jge	L443
	.p2align 4,,15
L499:
	mov	eax, DWORD PTR [ebp+16]
	mov	ecx, DWORD PTR [ebp+36]
	mov	edx, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+esi*8]
	mov	eax, DWORD PTR [ebp+32]
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [esp+8]
	fld	QWORD PTR [edx+esi*8]
	fsubr	QWORD PTR [eax+ebx*8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fld	QWORD PTR [ebp+52]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jae	L488
	fdiv	QWORD PTR [ebp+52]
	fld	st(0)
	fmul	DWORD PTR LC128
	fmulp	st(1), st
	fadd	QWORD PTR [ebp-32]
L395:
	mov	eax, DWORD PTR [ebp+20]
	mov	ecx, DWORD PTR [ebp+40]
	mov	edx, DWORD PTR [ebp+24]
	fld	QWORD PTR [eax+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [edx+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [ecx+ebx*8]
	fld	QWORD PTR [edi+ebx*8]
	fsubrp	st(1), st
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jl	L499
L443:
	inc	esi
	jmp	L473
	.p2align 4,,7
L488:
	fstp	QWORD PTR [esp]
	call	_log
	jmp	L395
L481:
	xor	esi, esi
	.p2align 4,,15
L471:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L373
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jge	L439
	.p2align 4,,15
L500:
	mov	edx, DWORD PTR [ebp+32]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [edx+ebx*8]
	fld	QWORD PTR [eax+esi*8]
	fxch	st(1)
	fucomi	st, st(1)
	jne	L451
	jp	L451
	mov	eax, DWORD PTR [ebp+36]
	mov	ecx, DWORD PTR [ebp+16]
	fld	QWORD PTR [eax+ebx*8]
	fld	QWORD PTR [ecx+esi*8]
	fxch	st(1)
	fucomi	st, st(1)
	jp	L491
	je	L492
L491:
	fxch	st(1)
	jmp	L383
	.p2align 4,,7
L492:
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	inc	ebx
	.p2align 4,,15
L501:
	cmp	ebx, DWORD PTR [ebp+28]
	jl	L500
L439:
	inc	esi
	jmp	L471
	.p2align 4,,7
L451:
	mov	ecx, DWORD PTR [ebp+36]
	mov	edx, DWORD PTR [ebp+16]
	fld	QWORD PTR [ecx+ebx*8]
	fld	QWORD PTR [edx+esi*8]
L383:
	fsubp	st(1), st
	fxch	st(1)
	fsubrp	st(2), st
	fstp	QWORD PTR [esp+8]
	fstp	QWORD PTR [esp]
	call	_hypot
	fstp	QWORD PTR [esp]
	call	_log
	mov	ecx, DWORD PTR [ebp+20]
	mov	eax, DWORD PTR [ebp+24]
	mov	edx, DWORD PTR [ebp+40]
	fld	QWORD PTR [ecx+esi*8]
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [eax+esi*8]
	fxch	st(1)
	fsubr	QWORD PTR [edx+ebx*8]
	fxch	st(1)
	fsubr	QWORD PTR [edi+ebx*8]
	fxch	st(1)
	fstp	QWORD PTR [edx+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	inc	ebx
	jmp	L501
L485:
	fst	QWORD PTR [esp]
	fld	st(0)
	fld	st(1)
	fxch	st(1)
	fmul	QWORD PTR LC131
	fxch	st(1)
	fmul	QWORD PTR LC132
	fxch	st(1)
	fadd	QWORD PTR LC133
	fxch	st(1)
	fadd	QWORD PTR LC134
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC135
	fxch	st(1)
	fadd	QWORD PTR LC136
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC137
	fxch	st(1)
	fadd	QWORD PTR LC138
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fadd	QWORD PTR LC139
	fxch	st(1)
	fadd	QWORD PTR LC140
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fsub	QWORD PTR LC141
	fxch	st(1)
	fsub	QWORD PTR LC142
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fadd	QWORD PTR LC143
	fxch	st(1)
	fsub	QWORD PTR LC144
	fdivp	st(1), st
	fstp	QWORD PTR [ebp-104]
	call	_log
	fld	QWORD PTR [ebp-104]
	fsubrp	st(1), st
	jmp	L415
L486:
	fld	QWORD PTR [ebp+60]
	fstp	QWORD PTR [esp]
	call	_log
	fadd	st, st(0)
	fstp	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp+60]
	fmul	st, st(0)
	fadd	QWORD PTR [ebp-72]
	fstp	QWORD PTR [esp]
	call	_log
	fdivr	QWORD PTR [ebp-40]
	fstp	QWORD PTR [ebp-40]
	jmp	L399
	.section .rdata,"dr"
	.align 8
LC171:
	.long	-1670671864
	.long	1072962135
	.align 4
LC172:
	.long	-1082130432
	.text
	.align 2
	.p2align 4,,15
.globl __Z4zinviPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb
	.def	__Z4zinviPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z4zinviPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 972
	mov	eax, DWORD PTR [ebp+48]
	fld	QWORD PTR [ebp+52]
	fxch	st(1)
	mov	edi, DWORD PTR [ebp+40]
	cmp	eax, 1
	movzx	edx, BYTE PTR [ebp+68]
	fstp	QWORD PTR [ebp-736]
	je	L528
	jle	L701
	cmp	eax, 2
	je	L615
	fstp	st(0)
	cmp	eax, 3
	je	L702
L730:
	add	esp, 972
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L528:
	fld	st(0)
	xor	ecx, ecx
	fmul	st, st(1)
	cmp	ecx, DWORD PTR [ebp+8]
	fdivr	QWORD PTR [ebp-736]
	fstp	QWORD PTR [ebp-744]
	jge	L731
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	ebx, [ebp-688]
	fld	QWORD PTR _I+8
	fxch	st(1)
	lea	esi, [ebp-696]
	fstp	QWORD PTR [ebp-864]
	fstp	QWORD PTR [ebp-872]
	.p2align 4,,15
L569:
	xor	edx, edx
	cmp	edx, DWORD PTR [ebp+28]
	jge	L669
	movsd	QWORD PTR [ebp-960], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-960]
	fld	QWORD PTR [ebp-864]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+16]
	fmul	st, st(2)
	fxch	st(1)
	fstp	QWORD PTR [ebp-816]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-752]
	fstp	QWORD PTR [ebp-824]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+24]
	fstp	QWORD PTR [ebp-832]
	fld	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [ebp-840]
	jmp	L568
	.p2align 4,,7
L712:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-200]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-192]
	jae	L703
L732:
	fld	QWORD PTR [ebp-128]
	fld	QWORD PTR [ebp-136]
	fld	QWORD PTR [ebp-744]
	fxch	st(2)
	fchs
	fxch	st(2)
	fmul	st, st(1)
	fld	st(2)
	fst	QWORD PTR [ebp-272]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fst	QWORD PTR [ebp-256]
	fmul	QWORD PTR [ebp-744]
	fxch	st(1)
	fsubrp	st(3), st
	fxch	st(1)
	fst	QWORD PTR [ebp-280]
	fst	QWORD PTR [ebp-264]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-296]
	fst	QWORD PTR [ebp-248]
	fstp	QWORD PTR [ebp-200]
	faddp	st(1), st
	fst	QWORD PTR [ebp-288]
	fst	QWORD PTR [ebp-240]
L727:
	fstp	QWORD PTR [ebp-192]
	fld	QWORD PTR [ebp-200]
	mov	eax, DWORD PTR [ebp+44]
	fld	QWORD PTR [ebp-192]
	fld	QWORD PTR [ebp-832]
	fld	QWORD PTR [ebp-840]
	fxch	st(1)
	fmul	st, st(3)
	fxch	st(1)
	fmul	st, st(2)
	fxch	st(3)
	fmul	QWORD PTR [ebp-840]
	fxch	st(2)
	fmul	QWORD PTR [ebp-832]
	fxch	st(1)
	fsubrp	st(3), st
	fxch	st(2)
	fsubr	QWORD PTR [edi+edx*8]
	fxch	st(2)
	faddp	st(1), st
	fsubr	QWORD PTR [eax+edx*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+28]
	jge	L711
L568:
	fld	QWORD PTR [ebp-816]
	mov	eax, DWORD PTR [ebp+32]
	fld	QWORD PTR [ebp-824]
	fld	QWORD PTR [ebp-864]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+36]
	fld	QWORD PTR [ebp-872]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(3)
	mov	eax, esi
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-872]
	fxch	st(2)
	fsubrp	st(3), st
	fxch	st(1)
	fadd	QWORD PTR [ebp-752]
	fxch	st(1)
	fadd	st, st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-168]
	fstp	QWORD PTR [ebp-152]
	fld	st(0)
	fld	st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-160]
	fxch	st(2)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fst	QWORD PTR [ebp-144]
	fst	QWORD PTR [ebp-176]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	QWORD PTR [ebp-696]
	fxch	st(2)
	cmovbe	eax, ebx
	fst	QWORD PTR [ebp-184]
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fst	QWORD PTR [ebp-128]
	fxch	st(2)
	fstp	QWORD PTR [ebp-688]
	fld	QWORD PTR [eax]
	fucomi	st, st(3)
	fld	st(0)
	jp	L712
	jne	L712
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fxch	st(1)
	fst	QWORD PTR [ebp-200]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-192]
	jb	L732
L703:
	fld	QWORD PTR [ebp-136]
	fld	QWORD PTR [ebp-128]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fld	st(4)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-232]
	fxch	st(1)
	fst	QWORD PTR [ebp-224]
	fxch	st(1)
	fst	QWORD PTR [ebp-216]
	fxch	st(1)
	fst	QWORD PTR [ebp-208]
	fxch	st(1)
	fstp	QWORD PTR [ebp-200]
	jmp	L727
L615:
	fmul	st, st(0)
	test	dl, dl
	fdivr	QWORD PTR LC171
	fst	QWORD PTR [ebp-776]
	jne	L704
	fstp	st(0)
L616:
	xor	esi, esi
	.p2align 4,,15
L699:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L730
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jge	L676
	fldz
	jmp	L659
	.p2align 4,,7
L715:
	fstp	st(0)
	fdiv	st(2), st
	fdiv	st(1), st
	fxch	st(2)
	fmul	st, st(0)
	fxch	st(1)
	fmul	st, st(0)
	faddp	st(1), st
	fsqrt
	fmulp	st(1), st
L636:
	fucomi	st, st(1)
	fld	st(1)
	jp	L698
	je	L717
L698:
	fld	QWORD PTR [ebp+60]
	fxch	st(1)
	fst	QWORD PTR [ebp-584]
	fst	QWORD PTR [ebp-576]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L733
	fstp	st(0)
	fld	QWORD PTR [ebp-520]
	fld	QWORD PTR [ebp-512]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fxch	st(4)
	fmul	st, st(3)
	fxch	st(4)
	faddp	st(3), st
	fxch	st(3)
	fsubrp	st(1), st
	fxch	st(1)
	fdiv	st, st(2)
	fxch	st(1)
	fdivrp	st(2), st
	fst	QWORD PTR [ebp-616]
	fxch	st(1)
	fst	QWORD PTR [ebp-608]
	fxch	st(1)
	fst	QWORD PTR [ebp-600]
	fxch	st(1)
	fst	QWORD PTR [ebp-592]
	fxch	st(1)
	fstp	QWORD PTR [ebp-584]
L728:
	fstp	QWORD PTR [ebp-576]
	fld	QWORD PTR [ebp-584]
	mov	eax, DWORD PTR [ebp+20]
	mov	edx, DWORD PTR [ebp+24]
	fld	QWORD PTR [ebp-576]
	mov	ecx, DWORD PTR [ebp+44]
	fld	QWORD PTR [eax+esi*8]
	fld	QWORD PTR [edx+esi*8]
	fld	st(1)
	fld	st(1)
	fxch	st(1)
	fmul	st, st(5)
	fxch	st(1)
	fmul	st, st(4)
	fxch	st(3)
	fmulp	st(4), st
	fxch	st(1)
	fmulp	st(4), st
	fsubrp	st(1), st
	fsubr	QWORD PTR [edi+ebx*8]
	fxch	st(1)
	faddp	st(2), st
	fxch	st(1)
	fsubr	QWORD PTR [ecx+ebx*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [ecx+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jge	L714
L659:
	fld	QWORD PTR _I
	mov	eax, DWORD PTR [ebp+12]
	mov	ecx, DWORD PTR [ebp+16]
	fld	QWORD PTR _I+8
	mov	edx, DWORD PTR [ebp+32]
	fld	st(1)
	fld	QWORD PTR [eax+esi*8]
	fxch	st(3)
	mov	eax, DWORD PTR [ebp+36]
	fmul	st, st(4)
	fld	QWORD PTR [ecx+esi*8]
	lea	ecx, [ebp-728]
	fld	QWORD PTR [edx+ebx*8]
	fxch	st(1)
	lea	edx, [ebp-720]
	fsubr	QWORD PTR [eax+ebx*8]
	fxch	st(1)
	fsubrp	st(5), st
	fld	st(3)
	fmul	st, st(6)
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(1)
	fmulp	st(4), st
	fsubrp	st(2), st
	faddp	st(2), st
	fadd	st(2), st
	fst	QWORD PTR [ebp-552]
	fld	st(1)
	fld	st(3)
	fxch	st(3)
	fst	QWORD PTR [ebp-544]
	fxch	st(3)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fstp	QWORD PTR [ebp-536]
	fxch	st(2)
	fst	QWORD PTR [ebp-528]
	fxch	st(1)
	fucomi	st, st(2)
	fxch	st(2)
	fstp	QWORD PTR [ebp-720]
	cmovbe	ecx, edx
	fst	QWORD PTR [ebp-560]
	fxch	st(2)
	fst	QWORD PTR [ebp-568]
	fst	QWORD PTR [ebp-520]
	fxch	st(2)
	fst	QWORD PTR [ebp-512]
	fxch	st(1)
	fstp	QWORD PTR [ebp-728]
	fld	QWORD PTR [ecx]
	fucomi	st, st(3)
	fld	st(0)
	jp	L715
	jne	L715
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	jmp	L636
	.p2align 4,,7
L711:
	fstp	st(0)
L669:
	inc	ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jl	L569
L731:
	fstp	st(0)
	add	esp, 972
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 4,,7
L717:
	fstp	st(0)
	fstp	st(0)
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+28]
	jl	L659
L714:
	fstp	st(0)
L676:
	inc	esi
	jmp	L699
L733:
	fstp	st(1)
	fld	QWORD PTR [ebp-776]
	fchs
	fmul	st, st(1)
	fmulp	st(1), st
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-952]
	call	_expm1
	fld	QWORD PTR [ebp-520]
	fld	QWORD PTR [ebp-512]
	fxch	st(2)
	fchs
	fld	st(1)
	fld	st(3)
	fld	QWORD PTR [ebp-952]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(1)
	fmul	st, st(5)
	fld	st(3)
	fmul	st, st(5)
	fxch	st(2)
	faddp	st(1), st
	fxch	st(3)
	fmul	st, st(5)
	fxch	st(4)
	fmul	st, st(2)
	fxch	st(5)
	fmul	st, st(2)
	fxch	st(5)
	fsubrp	st(4), st
	faddp	st(4), st
	fld	QWORD PTR [ebp-736]
	fxch	st(3)
	fdiv	st, st(2)
	fxch	st(4)
	fdivrp	st(2), st
	fxch	st(3)
	fst	QWORD PTR [ebp-656]
	fld	st(0)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(1)
	fst	QWORD PTR [ebp-640]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fmul	QWORD PTR [ebp-736]
	fxch	st(2)
	fst	QWORD PTR [ebp-664]
	fxch	st(1)
	fsubrp	st(3), st
	fst	QWORD PTR [ebp-648]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-680]
	fst	QWORD PTR [ebp-632]
	fstp	QWORD PTR [ebp-584]
	faddp	st(1), st
	fst	QWORD PTR [ebp-672]
	fst	QWORD PTR [ebp-624]
	jmp	L728
L701:
	test	eax, eax
	jne	L731
	fmul	st, st(0)
	test	dl, dl
	fstp	QWORD PTR [ebp-760]
	jne	L706
L571:
	xor	ecx, ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jge	L730
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	ebx, [ebp-704]
	fld	QWORD PTR _I+8
	fxch	st(1)
	lea	esi, [ebp-712]
	fstp	QWORD PTR [ebp-848]
	fstp	QWORD PTR [ebp-856]
L614:
	xor	edx, edx
	cmp	edx, DWORD PTR [ebp+28]
	jge	L673
	movsd	QWORD PTR [ebp-960], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-960]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+16]
	fstp	QWORD PTR [ebp-792]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-800]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+24]
	fstp	QWORD PTR [ebp-808]
	fld	QWORD PTR [ebp-848]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	fmul	st, st(2)
	fstp	QWORD PTR [ebp-768]
	jmp	L613
	.p2align 4,,7
L720:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fld	QWORD PTR [ebp+60]
	fxch	st(3)
	fst	QWORD PTR [ebp-376]
	fst	QWORD PTR [ebp-368]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jae	L707
L734:
	fld	QWORD PTR [ebp-304]
	fxch	st(3)
	fmul	st, st(0)
	fld	QWORD PTR [ebp-312]
	fxch	st(4)
	fchs
	fxch	st(1)
	fadd	QWORD PTR [ebp-760]
	fld	st(4)
	fld	st(2)
	fxch	st(6)
	fst	QWORD PTR [ebp-472]
	fld	st(2)
	fmul	st, st(3)
	fxch	st(4)
	fst	QWORD PTR [ebp-464]
	fxch	st(2)
	fmul	st, st(3)
	fxch	st(7)
	fmul	st, st(5)
	fxch	st(1)
	fst	QWORD PTR [ebp-456]
	fmul	st, st(5)
	fxch	st(4)
	fadd	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-448]
	fmulp	st(3), st
	faddp	st(6), st
	fld	QWORD PTR [ebp-736]
	fxch	st(2)
	fsubrp	st(3), st
	fdiv	st(5), st
	fdivp	st(2), st
	fld	st(1)
	fxch	st(1)
	fmul	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-480]
	fxch	st(1)
	fmul	st, st(3)
	fxch	st(1)
	fst	QWORD PTR [ebp-432]
	fmul	QWORD PTR [ebp-736]
	fxch	st(2)
	fsubrp	st(1), st
	fxch	st(4)
	fst	QWORD PTR [ebp-488]
	fst	QWORD PTR [ebp-440]
	fmul	st, st(2)
	fxch	st(4)
	fst	QWORD PTR [ebp-504]
	fst	QWORD PTR [ebp-424]
	fstp	QWORD PTR [ebp-376]
	faddp	st(3), st
	fxch	st(2)
	fst	QWORD PTR [ebp-496]
	fst	QWORD PTR [ebp-416]
L729:
	fstp	QWORD PTR [ebp-368]
	fld	QWORD PTR [ebp-368]
	fld	st(1)
	mov	eax, DWORD PTR [ebp+44]
	fld	QWORD PTR [ebp-376]
	fxch	st(1)
	fmul	st, st(2)
	fld	QWORD PTR [ebp-808]
	fxch	st(3)
	fmul	QWORD PTR [ebp-808]
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(2)
	fsubrp	st(1), st
	fxch	st(2)
	faddp	st(1), st
	fxch	st(1)
	fsubr	QWORD PTR [edi+edx*8]
	fxch	st(1)
	fsubr	QWORD PTR [eax+edx*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+28]
	jge	L719
L613:
	fld	QWORD PTR [ebp-792]
	mov	eax, DWORD PTR [ebp+32]
	fld	QWORD PTR [ebp-800]
	fld	QWORD PTR [ebp-848]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+36]
	fld	QWORD PTR [ebp-856]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(5)
	fxch	st(3)
	mov	eax, esi
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-856]
	fxch	st(2)
	fsubrp	st(3), st
	fxch	st(1)
	fadd	QWORD PTR [ebp-768]
	fxch	st(1)
	fadd	st, st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-344]
	fstp	QWORD PTR [ebp-328]
	fld	st(0)
	fld	st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-336]
	fxch	st(2)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fst	QWORD PTR [ebp-320]
	fst	QWORD PTR [ebp-352]
	fxch	st(2)
	fucomi	st, st(1)
	fxch	st(1)
	fstp	QWORD PTR [ebp-704]
	fxch	st(2)
	cmovbe	eax, ebx
	fst	QWORD PTR [ebp-360]
	fst	QWORD PTR [ebp-312]
	fxch	st(1)
	fst	QWORD PTR [ebp-304]
	fxch	st(2)
	fstp	QWORD PTR [ebp-712]
	fld	QWORD PTR [eax]
	fucomi	st, st(4)
	fld	st(0)
	jp	L720
	jne	L720
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fld	QWORD PTR [ebp+60]
	fxch	st(3)
	fst	QWORD PTR [ebp-376]
	fst	QWORD PTR [ebp-368]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jb	L734
L707:
	fstp	st(2)
	fld	QWORD PTR [ebp-312]
	fld	QWORD PTR [ebp-304]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(6)
	fxch	st(1)
	faddp	st(2), st
	fld	st(5)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-408]
	fxch	st(1)
	fst	QWORD PTR [ebp-400]
	fxch	st(1)
	fst	QWORD PTR [ebp-392]
	fxch	st(1)
	fst	QWORD PTR [ebp-384]
	fxch	st(1)
	fstp	QWORD PTR [ebp-376]
	jmp	L729
L702:
	xor	edx, edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L730
	fld	QWORD PTR [ebp-120]
	fld	QWORD PTR _I
	fld	QWORD PTR _I+8
	fld	QWORD PTR [ebp-80]
	fxch	st(3)
	fstp	QWORD PTR [ebp-904]
	fld	QWORD PTR [ebp-112]
	fxch	st(2)
	fstp	QWORD PTR [ebp-888]
	fstp	QWORD PTR [ebp-896]
	fld	QWORD PTR [ebp-72]
	fxch	st(1)
	fstp	QWORD PTR [ebp-912]
	fld	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp-88]
	fld	QWORD PTR [ebp-64]
	fld	QWORD PTR [ebp-104]
	fxch	st(3)
	fstp	QWORD PTR [ebp-920]
	fld	QWORD PTR [ebp-32]
	fld	QWORD PTR [ebp-96]
	fxch	st(1)
	movsd	xmm0, QWORD PTR [ebp-56]
	movsd	xmm1, QWORD PTR [ebp-48]
	fstp	QWORD PTR [ebp-928]
	xor	eax, eax
	cmp	eax, DWORD PTR [ebp+28]
	jge	L722
L737:
	fxch	st(2)
	fstp	QWORD PTR [ebp-880]
	mov	ebx, DWORD PTR [ebp+12]
	movapd	xmm4, xmm0
	movsd	xmm0, QWORD PTR [ebp-928]
	movapd	xmm3, xmm1
	mov	ecx, DWORD PTR [ebp+16]
	movsd	xmm2, QWORD PTR [ebx+edx*8]
	movsd	xmm1, QWORD PTR [ebp-920]
	fld	QWORD PTR [ecx+edx*8]
	fxch	st(4)
	fstp	QWORD PTR [ebp-960]
	movsd	xmm7, QWORD PTR [ebp-960]
	fstp	QWORD PTR [ebp-960]
	movsd	xmm6, QWORD PTR [ebp-960]
	fstp	QWORD PTR [ebp-960]
	movsd	xmm5, QWORD PTR [ebp-960]
	fld	QWORD PTR [ebp-904]
	fld	QWORD PTR [ebp-912]
	fxch	st(3)
	fstp	QWORD PTR [ebp-784]
	jmp	L526
	.p2align 4,,7
L723:
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	mov	ebx, DWORD PTR [ebp+36]
	fld	QWORD PTR [ebx+eax*8]
L512:
	movsd	QWORD PTR [ebp-960], xmm2
	fsub	QWORD PTR [ebp-784]
	fldz
	mov	ebx, DWORD PTR [ebp+20]
	fld	QWORD PTR [ebp-960]
	fldz
	fxch	st(2)
	mov	ecx, DWORD PTR [ebp+24]
	fmul	QWORD PTR [ebp-888]
	fxch	st(1)
	mov	esi, DWORD PTR [ebp+44]
	fsubp	st(4), st
	fld	QWORD PTR [ebp-888]
	fldz
	fxch	st(1)
	fmul	st, st(4)
	fxch	st(4)
	fmul	QWORD PTR [ebp-896]
	fxch	st(4)
	fstp	QWORD PTR [ebp-880]
	fldz
	fxch	st(2)
	faddp	st(4), st
	fxch	st(1)
	fmul	QWORD PTR [ebp-896]
	fld	st(3)
	fmul	st, st(4)
	fxch	st(1)
	fsubr	QWORD PTR [ebp-880]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(4)
	fst	QWORD PTR [ebp-960]
	fxch	st(5)
	fadd	st, st(3)
	fxch	st(3)
	fstp	QWORD PTR [ebp-880]
	movsd	xmm6, QWORD PTR [ebp-960]
	fld	st(2)
	movsd	xmm7, QWORD PTR [ebp-880]
	fmul	st, st(3)
	fxch	st(2)
	fmul	st, st(3)
	fxch	st(4)
	movapd	xmm5, xmm6
	movapd	xmm3, xmm6
	fadd	st, st(3)
	fxch	st(3)
	fst	QWORD PTR [ebp-960]
	fxch	st(2)
	faddp	st(1), st
	fxch	st(3)
	fsub	st, st(4)
	fxch	st(2)
	movsd	xmm4, QWORD PTR [ebp-960]
	fdiv	st, st(3)
	fxch	st(2)
	fdivrp	st(3), st
	fxch	st(1)
	fst	QWORD PTR [ebp-960]
	fxch	st(2)
	movsd	xmm1, QWORD PTR [ebp-960]
	fst	QWORD PTR [ebp-960]
	fld	QWORD PTR [ebx+edx*8]
	fld	QWORD PTR [ecx+edx*8]
	movsd	xmm0, QWORD PTR [ebp-960]
	fld	st(1)
	fld	st(1)
	fxch	st(1)
	fmul	st, st(6)
	fxch	st(1)
	fmul	st, st(4)
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(2)
	fmul	st, st(6)
	fxch	st(1)
	fsubrp	st(3), st
	fxch	st(2)
	fsubr	QWORD PTR [edi+eax*8]
	fxch	st(1)
	faddp	st(2), st
	fxch	st(1)
	fsubr	QWORD PTR [esi+eax*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+eax*8]
	fstp	QWORD PTR [esi+eax*8]
	inc	eax
	cmp	eax, DWORD PTR [ebp+28]
	jge	L708
L735:
	fxch	st(2)
L526:
	mov	esi, DWORD PTR [ebp+32]
	fld	QWORD PTR [esi+eax*8]
	movsd	QWORD PTR [ebp-960], xmm2
	fld	QWORD PTR [ebp-960]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jne	L723
	jp	L723
	fld	QWORD PTR [ebp-784]
	mov	ecx, DWORD PTR [ebp+36]
	fld	QWORD PTR [ecx+eax*8]
	fucomi	st, st(1)
	fstp	st(1)
	jp	L724
	je	L725
L724:
	fstp	st(2)
	fstp	st(2)
	fstp	st(2)
	fstp	st(2)
	jmp	L512
L725:
	fstp	st(0)
	fstp	st(0)
	fxch	st(2)
	inc	eax
	cmp	eax, DWORD PTR [ebp+28]
	jl	L735
L708:
	movsd	QWORD PTR [ebp-960], xmm5
	fld	QWORD PTR [ebp-880]
	fxch	st(3)
	movsd	QWORD PTR [ebp-928], xmm0
	movsd	QWORD PTR [ebp-920], xmm1
	movapd	xmm0, xmm4
	fstp	QWORD PTR [ebp-904]
	movapd	xmm1, xmm3
	fld	QWORD PTR [ebp-960]
	fxch	st(1)
	movsd	QWORD PTR [ebp-960], xmm6
	fstp	QWORD PTR [ebp-912]
	fld	QWORD PTR [ebp-960]
	movsd	QWORD PTR [ebp-960], xmm7
	fld	QWORD PTR [ebp-960]
L665:
	inc	edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L736
	fxch	st(4)
	fxch	st(2)
	xor	eax, eax
	cmp	eax, DWORD PTR [ebp+28]
	jl	L737
L722:
	fxch	st(2)
	fxch	st(4)
	jmp	L665
L736:
	fld	QWORD PTR [ebp-928]
	fxch	st(3)
	movsd	QWORD PTR [ebp-48], xmm1
	movsd	QWORD PTR [ebp-56], xmm0
	fstp	QWORD PTR [ebp-96]
	fxch	st(2)
	fstp	QWORD PTR [ebp-32]
	fld	QWORD PTR [ebp-920]
	fxch	st(3)
	fstp	QWORD PTR [ebp-104]
	fstp	QWORD PTR [ebp-64]
	fstp	QWORD PTR [ebp-72]
	fstp	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp-912]
	fxch	st(2)
	fstp	QWORD PTR [ebp-80]
	fstp	QWORD PTR [ebp-88]
	fstp	QWORD PTR [ebp-112]
	fld	QWORD PTR [ebp-904]
	fstp	QWORD PTR [ebp-120]
	jmp	L730
L704:
	fchs
	fmul	QWORD PTR [ebp+60]
	fmul	QWORD PTR [ebp+60]
	fstp	QWORD PTR [esp]
	call	_expm1
	fld	DWORD PTR LC172
	fdivrp	st(1), st
	fstp	QWORD PTR [ebp-736]
	jmp	L616
L706:
	fld	QWORD PTR [ebp+60]
	fmul	st, st(0)
	fdivr	QWORD PTR [ebp-760]
	fadd	QWORD PTR [ebp-736]
	fstp	QWORD PTR [ebp-736]
	jmp	L571
L719:
	fstp	st(0)
	fstp	st(0)
L673:
	inc	ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jl	L614
	jmp	L730
	.section .rdata,"dr"
	.align 8
LC178:
	.long	-1670671864
	.long	1072962135
	.align 4
LC179:
	.long	-1082130432
	.text
	.align 2
	.p2align 4,,15
.globl __Z4rinviPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb
	.def	__Z4rinviPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z4rinviPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 956
	mov	eax, DWORD PTR [ebp+44]
	fld	QWORD PTR [ebp+48]
	fxch	st(1)
	mov	edi, DWORD PTR [ebp+36]
	cmp	eax, 1
	movzx	edx, BYTE PTR [ebp+64]
	fstp	QWORD PTR [ebp-736]
	je	L762
	jle	L929
	cmp	eax, 2
	je	L845
	fstp	st(0)
	cmp	eax, 3
	je	L930
L958:
	add	esp, 956
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L762:
	fld	st(0)
	xor	ecx, ecx
	fmul	st, st(1)
	cmp	ecx, DWORD PTR [ebp+8]
	fdivr	QWORD PTR [ebp-736]
	fstp	QWORD PTR [ebp-744]
	jge	L959
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	esi, [ebp-688]
	fld	QWORD PTR _I+8
	fxch	st(1)
	lea	ebx, [ebp-696]
	fstp	QWORD PTR [ebp-848]
	fstp	QWORD PTR [ebp-856]
	.p2align 4,,15
L801:
	xor	edx, edx
	cmp	edx, DWORD PTR [ebp+24]
	jge	L897
	movsd	QWORD PTR [ebp-928], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-928]
	fld	QWORD PTR [ebp-848]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+16]
	fmul	st, st(2)
	fxch	st(1)
	fstp	QWORD PTR [ebp-784]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-752]
	fstp	QWORD PTR [ebp-792]
	fld	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [ebp-800]
	jmp	L800
	.p2align 4,,7
L940:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-200]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-192]
	jae	L931
L960:
	fld	QWORD PTR [ebp-128]
	fld	QWORD PTR [ebp-136]
	fld	QWORD PTR [ebp-744]
	fxch	st(2)
	fchs
	fxch	st(2)
	fmul	st, st(1)
	fld	st(2)
	fst	QWORD PTR [ebp-272]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fst	QWORD PTR [ebp-256]
	fmul	QWORD PTR [ebp-744]
	fxch	st(1)
	fsubrp	st(3), st
	fxch	st(1)
	fst	QWORD PTR [ebp-280]
	fst	QWORD PTR [ebp-264]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-296]
	fst	QWORD PTR [ebp-248]
	fstp	QWORD PTR [ebp-200]
	faddp	st(1), st
	fst	QWORD PTR [ebp-288]
	fst	QWORD PTR [ebp-240]
L955:
	fstp	QWORD PTR [ebp-192]
	fld	QWORD PTR [ebp-800]
	mov	eax, DWORD PTR [ebp+40]
	fmul	QWORD PTR [ebp-200]
	fsubr	QWORD PTR [edi+edx*8]
	fstp	QWORD PTR [edi+edx*8]
	fld	QWORD PTR [ebp-800]
	fmul	QWORD PTR [ebp-192]
	fsubr	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+24]
	jge	L939
L800:
	fld	QWORD PTR [ebp-784]
	mov	eax, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp-792]
	fld	QWORD PTR [ebp-848]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+32]
	fld	QWORD PTR [ebp-856]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(3)
	mov	eax, ebx
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-856]
	fxch	st(2)
	fsubrp	st(3), st
	fxch	st(1)
	fadd	QWORD PTR [ebp-752]
	fxch	st(1)
	fadd	st, st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-168]
	fstp	QWORD PTR [ebp-152]
	fld	st(0)
	fld	st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-160]
	fxch	st(2)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fst	QWORD PTR [ebp-144]
	fst	QWORD PTR [ebp-176]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	QWORD PTR [ebp-696]
	fxch	st(2)
	cmovbe	eax, esi
	fst	QWORD PTR [ebp-184]
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fst	QWORD PTR [ebp-128]
	fxch	st(2)
	fstp	QWORD PTR [ebp-688]
	fld	QWORD PTR [eax]
	fucomi	st, st(3)
	fld	st(0)
	jp	L940
	jne	L940
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fxch	st(1)
	fst	QWORD PTR [ebp-200]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-192]
	jb	L960
L931:
	fld	QWORD PTR [ebp-136]
	fld	QWORD PTR [ebp-128]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fld	st(4)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-232]
	fxch	st(1)
	fst	QWORD PTR [ebp-224]
	fxch	st(1)
	fst	QWORD PTR [ebp-216]
	fxch	st(1)
	fst	QWORD PTR [ebp-208]
	fxch	st(1)
	fstp	QWORD PTR [ebp-200]
	jmp	L955
L845:
	fmul	st, st(0)
	test	dl, dl
	fdivr	QWORD PTR LC178
	fst	QWORD PTR [ebp-776]
	jne	L932
	fstp	st(0)
L846:
	xor	esi, esi
	.p2align 4,,15
L927:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L958
	xor	ebx, ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jge	L904
	fldz
	jmp	L887
	.p2align 4,,7
L943:
	fstp	st(0)
	fdiv	st(2), st
	fdiv	st(1), st
	fxch	st(2)
	fmul	st, st(0)
	fxch	st(1)
	fmul	st, st(0)
	faddp	st(1), st
	fsqrt
	fmulp	st(1), st
L866:
	fucomi	st, st(1)
	fld	st(1)
	jp	L926
	je	L945
L926:
	fld	QWORD PTR [ebp+56]
	fxch	st(1)
	fst	QWORD PTR [ebp-584]
	fst	QWORD PTR [ebp-576]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L961
	fstp	st(0)
	fld	QWORD PTR [ebp-520]
	fld	QWORD PTR [ebp-512]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fxch	st(4)
	fmul	st, st(3)
	fxch	st(4)
	faddp	st(3), st
	fxch	st(3)
	fsubrp	st(1), st
	fxch	st(1)
	fdiv	st, st(2)
	fxch	st(1)
	fdivrp	st(2), st
	fst	QWORD PTR [ebp-616]
	fxch	st(1)
	fst	QWORD PTR [ebp-608]
	fxch	st(1)
	fst	QWORD PTR [ebp-600]
	fxch	st(1)
	fst	QWORD PTR [ebp-592]
	fxch	st(1)
	fstp	QWORD PTR [ebp-584]
L956:
	fstp	QWORD PTR [ebp-576]
	mov	edx, DWORD PTR [ebp+20]
	mov	ecx, DWORD PTR [ebp+40]
	fld	QWORD PTR [edx+esi*8]
	fld	st(0)
	fmul	QWORD PTR [ebp-584]
	fsubr	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [edi+ebx*8]
	fmul	QWORD PTR [ebp-576]
	fsubr	QWORD PTR [ecx+ebx*8]
	fstp	QWORD PTR [ecx+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jge	L942
L887:
	fld	QWORD PTR _I
	mov	eax, DWORD PTR [ebp+12]
	mov	ecx, DWORD PTR [ebp+16]
	fld	QWORD PTR _I+8
	mov	edx, DWORD PTR [ebp+28]
	fld	st(1)
	fld	QWORD PTR [eax+esi*8]
	fxch	st(3)
	mov	eax, DWORD PTR [ebp+32]
	fmul	st, st(4)
	fld	QWORD PTR [ecx+esi*8]
	lea	ecx, [ebp-728]
	fld	QWORD PTR [edx+ebx*8]
	fxch	st(1)
	lea	edx, [ebp-720]
	fsubr	QWORD PTR [eax+ebx*8]
	fxch	st(1)
	fsubrp	st(5), st
	fld	st(3)
	fmul	st, st(6)
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(1)
	fmulp	st(4), st
	fsubrp	st(2), st
	faddp	st(2), st
	fadd	st(2), st
	fst	QWORD PTR [ebp-552]
	fld	st(1)
	fld	st(3)
	fxch	st(3)
	fst	QWORD PTR [ebp-544]
	fxch	st(3)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fstp	QWORD PTR [ebp-536]
	fxch	st(2)
	fst	QWORD PTR [ebp-528]
	fxch	st(1)
	fucomi	st, st(2)
	fxch	st(2)
	fstp	QWORD PTR [ebp-720]
	cmovbe	ecx, edx
	fst	QWORD PTR [ebp-560]
	fxch	st(2)
	fst	QWORD PTR [ebp-568]
	fst	QWORD PTR [ebp-520]
	fxch	st(2)
	fst	QWORD PTR [ebp-512]
	fxch	st(1)
	fstp	QWORD PTR [ebp-728]
	fld	QWORD PTR [ecx]
	fucomi	st, st(3)
	fld	st(0)
	jp	L943
	jne	L943
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	jmp	L866
	.p2align 4,,7
L939:
	fstp	st(0)
L897:
	inc	ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jl	L801
L959:
	fstp	st(0)
	add	esp, 956
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 4,,7
L945:
	fstp	st(0)
	fstp	st(0)
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+24]
	jl	L887
L942:
	fstp	st(0)
L904:
	inc	esi
	jmp	L927
L961:
	fstp	st(1)
	fld	QWORD PTR [ebp-776]
	fchs
	fmul	st, st(1)
	fmulp	st(1), st
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-920]
	call	_expm1
	fld	QWORD PTR [ebp-520]
	fld	QWORD PTR [ebp-512]
	fxch	st(2)
	fchs
	fld	st(1)
	fld	st(3)
	fld	QWORD PTR [ebp-920]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(1)
	fmul	st, st(5)
	fld	st(3)
	fmul	st, st(5)
	fxch	st(2)
	faddp	st(1), st
	fxch	st(3)
	fmul	st, st(5)
	fxch	st(4)
	fmul	st, st(2)
	fxch	st(5)
	fmul	st, st(2)
	fxch	st(5)
	fsubrp	st(4), st
	faddp	st(4), st
	fld	QWORD PTR [ebp-736]
	fxch	st(3)
	fdiv	st, st(2)
	fxch	st(4)
	fdivrp	st(2), st
	fxch	st(3)
	fst	QWORD PTR [ebp-656]
	fld	st(0)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(1)
	fst	QWORD PTR [ebp-640]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fmul	QWORD PTR [ebp-736]
	fxch	st(2)
	fst	QWORD PTR [ebp-664]
	fxch	st(1)
	fsubrp	st(3), st
	fst	QWORD PTR [ebp-648]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-680]
	fst	QWORD PTR [ebp-632]
	fstp	QWORD PTR [ebp-584]
	faddp	st(1), st
	fst	QWORD PTR [ebp-672]
	fst	QWORD PTR [ebp-624]
	jmp	L956
L929:
	test	eax, eax
	jne	L959
	fmul	st, st(0)
	test	dl, dl
	fstp	QWORD PTR [ebp-760]
	jne	L934
L803:
	xor	ecx, ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jge	L958
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	ebx, [ebp-704]
	fld	QWORD PTR _I+8
	fxch	st(1)
	lea	esi, [ebp-712]
	fstp	QWORD PTR [ebp-832]
	fstp	QWORD PTR [ebp-840]
L844:
	xor	edx, edx
	cmp	edx, DWORD PTR [ebp+24]
	jge	L901
	movsd	QWORD PTR [ebp-928], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-928]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+16]
	fstp	QWORD PTR [ebp-808]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-816]
	fld	QWORD PTR [ebp-832]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	fmul	st, st(2)
	fstp	QWORD PTR [ebp-768]
	jmp	L843
	.p2align 4,,7
L948:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fld	QWORD PTR [ebp+56]
	fxch	st(3)
	fst	QWORD PTR [ebp-376]
	fst	QWORD PTR [ebp-368]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jae	L935
L962:
	fld	QWORD PTR [ebp-304]
	fxch	st(3)
	fmul	st, st(0)
	fld	QWORD PTR [ebp-312]
	fxch	st(4)
	fchs
	fxch	st(1)
	fadd	QWORD PTR [ebp-760]
	fld	st(4)
	fld	st(2)
	fxch	st(6)
	fst	QWORD PTR [ebp-472]
	fld	st(2)
	fmul	st, st(3)
	fxch	st(4)
	fst	QWORD PTR [ebp-464]
	fxch	st(2)
	fmul	st, st(3)
	fxch	st(7)
	fmul	st, st(5)
	fxch	st(1)
	fst	QWORD PTR [ebp-456]
	fmul	st, st(5)
	fxch	st(4)
	fadd	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-448]
	fmulp	st(3), st
	faddp	st(6), st
	fld	QWORD PTR [ebp-736]
	fxch	st(2)
	fsubrp	st(3), st
	fdiv	st(5), st
	fdivp	st(2), st
	fld	st(1)
	fxch	st(1)
	fmul	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-480]
	fxch	st(1)
	fmul	st, st(3)
	fxch	st(1)
	fst	QWORD PTR [ebp-432]
	fmul	QWORD PTR [ebp-736]
	fxch	st(2)
	fsubrp	st(1), st
	fxch	st(4)
	fst	QWORD PTR [ebp-488]
	fst	QWORD PTR [ebp-440]
	fmul	st, st(2)
	fxch	st(4)
	fst	QWORD PTR [ebp-504]
	fst	QWORD PTR [ebp-424]
	fstp	QWORD PTR [ebp-376]
	faddp	st(3), st
	fxch	st(2)
	fst	QWORD PTR [ebp-496]
	fst	QWORD PTR [ebp-416]
L957:
	fstp	QWORD PTR [ebp-368]
	fld	st(0)
	mov	eax, DWORD PTR [ebp+40]
	fmul	QWORD PTR [ebp-376]
	fsubr	QWORD PTR [edi+edx*8]
	fstp	QWORD PTR [edi+edx*8]
	fld	st(0)
	fmul	QWORD PTR [ebp-368]
	fsubr	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+24]
	jge	L947
L843:
	fld	QWORD PTR [ebp-808]
	mov	eax, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp-816]
	fld	QWORD PTR [ebp-832]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+32]
	fld	QWORD PTR [ebp-840]
	fxch	st(2)
	fsubr	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(5)
	fxch	st(3)
	mov	eax, esi
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-840]
	fxch	st(2)
	fsubrp	st(3), st
	fxch	st(1)
	fadd	QWORD PTR [ebp-768]
	fxch	st(1)
	fadd	st, st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-344]
	fstp	QWORD PTR [ebp-328]
	fld	st(0)
	fld	st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-336]
	fxch	st(2)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fst	QWORD PTR [ebp-320]
	fst	QWORD PTR [ebp-352]
	fxch	st(2)
	fucomi	st, st(1)
	fxch	st(1)
	fstp	QWORD PTR [ebp-704]
	fxch	st(2)
	cmovbe	eax, ebx
	fst	QWORD PTR [ebp-360]
	fst	QWORD PTR [ebp-312]
	fxch	st(1)
	fst	QWORD PTR [ebp-304]
	fxch	st(2)
	fstp	QWORD PTR [ebp-712]
	fld	QWORD PTR [eax]
	fucomi	st, st(4)
	fld	st(0)
	jp	L948
	jne	L948
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fld	QWORD PTR [ebp+56]
	fxch	st(3)
	fst	QWORD PTR [ebp-376]
	fst	QWORD PTR [ebp-368]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jb	L962
L935:
	fstp	st(2)
	fld	QWORD PTR [ebp-312]
	fld	QWORD PTR [ebp-304]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(6)
	fxch	st(1)
	faddp	st(2), st
	fld	st(5)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-408]
	fxch	st(1)
	fst	QWORD PTR [ebp-400]
	fxch	st(1)
	fst	QWORD PTR [ebp-392]
	fxch	st(1)
	fst	QWORD PTR [ebp-384]
	fxch	st(1)
	fstp	QWORD PTR [ebp-376]
	jmp	L957
L930:
	xor	edx, edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L958
	fld	QWORD PTR _I+8
	fld	QWORD PTR _I
	fld	QWORD PTR [ebp-88]
	fld	QWORD PTR [ebp-80]
	fxch	st(3)
	fstp	QWORD PTR [ebp-872]
	fld	QWORD PTR [ebp-120]
	fxch	st(2)
	fstp	QWORD PTR [ebp-864]
	fld	QWORD PTR [ebp-72]
	fld	QWORD PTR [ebp-64]
	fxch	st(3)
	fstp	QWORD PTR [ebp-880]
	fld	QWORD PTR [ebp-112]
	fld	QWORD PTR [ebp-104]
	fxch	st(1)
	movsd	xmm0, QWORD PTR [ebp-96]
	movsd	xmm1, QWORD PTR [ebp-56]
	fstp	QWORD PTR [ebp-888]
	fld	QWORD PTR [ebp-40]
	movsd	xmm3, QWORD PTR [ebp-48]
	fstp	QWORD PTR [ebp-896]
	fld	QWORD PTR [ebp-32]
	fstp	QWORD PTR [ebp-904]
	xor	eax, eax
	cmp	eax, DWORD PTR [ebp+24]
	jge	L950
L965:
	mov	ecx, DWORD PTR [ebp+16]
	movapd	xmm4, xmm1
	movsd	xmm1, QWORD PTR [ebp-896]
	fldz
	mov	ebx, DWORD PTR [ebp+12]
	fld	QWORD PTR [ecx+edx*8]
	movapd	xmm5, xmm0
	fld	QWORD PTR [ebp-888]
	fxch	st(4)
	movsd	xmm2, QWORD PTR [ebx+edx*8]
	fstp	QWORD PTR [ebp-928]
	fxch	st(5)
	movsd	xmm7, QWORD PTR [ebp-928]
	fstp	QWORD PTR [ebp-928]
	fxch	st(4)
	fstp	QWORD PTR [ebp-824]
	movsd	xmm6, QWORD PTR [ebp-928]
	fld	QWORD PTR [ebp-904]
	fld	QWORD PTR [ebp-880]
	fxch	st(1)
	fstp	QWORD PTR [ebp-944]
	jmp	L760
	.p2align 4,,7
L951:
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fstp	st(2)
	mov	ebx, DWORD PTR [ebp+32]
	fld	QWORD PTR [ebx+eax*8]
	fstp	QWORD PTR [ebp-936]
L748:
	movsd	QWORD PTR [ebp-928], xmm2
	fld	QWORD PTR [ebp-936]
	mov	ecx, DWORD PTR [ebp+20]
	mov	esi, DWORD PTR [ebp+40]
	fld	QWORD PTR [ebp-928]
	fxch	st(1)
	fsub	QWORD PTR [ebp-824]
	fld	QWORD PTR [ebp-864]
	fxch	st(2)
	fsubp	st(4), st
	fld	QWORD PTR [ebp-872]
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(2)
	fmul	st, st(3)
	fld	QWORD PTR [ebp-864]
	fxch	st(2)
	fmul	QWORD PTR [ebp-872]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(3)
	fsubrp	st(1), st
	fxch	st(2)
	faddp	st(1), st
	fxch	st(3)
	fadd	st, st(1)
	fxch	st(1)
	fst	QWORD PTR [ebp-928]
	fld	st(3)
	movsd	xmm7, QWORD PTR [ebp-928]
	fld	st(2)
	fld	st(5)
	fst	QWORD PTR [ebp-928]
	fmul	st(2), st
	fxch	st(1)
	movsd	xmm6, QWORD PTR [ebp-928]
	fmul	st, st(4)
	fxch	st(6)
	fmul	st, st(5)
	movapd	xmm5, xmm6
	movapd	xmm3, xmm6
	fld	st(4)
	fst	QWORD PTR [ebp-928]
	fxch	st(7)
	faddp	st(3), st
	movsd	xmm4, QWORD PTR [ebp-928]
	fadd	st, st(6)
	fxch	st(4)
	fmul	st, st(5)
	fxch	st(4)
	fdiv	st, st(2)
	fxch	st(4)
	fsub	st, st(1)
	fxch	st(4)
	fst	QWORD PTR [ebp-928]
	fxch	st(4)
	movsd	xmm1, QWORD PTR [ebp-928]
	fdivrp	st(2), st
	fld	QWORD PTR [ecx+edx*8]
	fxch	st(2)
	fst	QWORD PTR [ebp-944]
	fld	st(2)
	fmul	st, st(5)
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(3)
	fsubr	QWORD PTR [edi+eax*8]
	fxch	st(3)
	fsubr	QWORD PTR [esi+eax*8]
	fxch	st(3)
	fstp	QWORD PTR [edi+eax*8]
	fxch	st(2)
	fstp	QWORD PTR [esi+eax*8]
	inc	eax
	cmp	eax, DWORD PTR [ebp+24]
	jge	L936
L963:
	fxch	st(5)
	fxch	st(1)
	fxch	st(2)
	fxch	st(3)
L760:
	mov	esi, DWORD PTR [ebp+28]
	fld	QWORD PTR [esi+eax*8]
	movsd	QWORD PTR [ebp-928], xmm2
	fld	QWORD PTR [ebp-928]
	fxch	st(1)
	fucomi	st, st(1)
	fstp	st(1)
	jne	L951
	jp	L951
	mov	ecx, DWORD PTR [ebp+32]
	fld	QWORD PTR [ecx+eax*8]
	fstp	QWORD PTR [ebp-936]
	movsd	xmm0, QWORD PTR [ebp-936]
	ucomisd	xmm0, QWORD PTR [ebp-824]
	jp	L952
	je	L953
L952:
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fstp	st(2)
	jmp	L748
L953:
	fstp	st(0)
	fxch	st(3)
	fxch	st(2)
	fxch	st(1)
	fxch	st(5)
	inc	eax
	cmp	eax, DWORD PTR [ebp+24]
	jl	L963
L936:
	fstp	st(4)
	movsd	QWORD PTR [ebp-928], xmm6
	fld	QWORD PTR [ebp-944]
	movapd	xmm0, xmm5
	fld	QWORD PTR [ebp-928]
	fxch	st(4)
	movsd	QWORD PTR [ebp-928], xmm7
	fstp	QWORD PTR [ebp-880]
	fld	QWORD PTR [ebp-928]
	fxch	st(1)
	fstp	QWORD PTR [ebp-904]
	fxch	st(1)
	movsd	QWORD PTR [ebp-896], xmm1
	movapd	xmm1, xmm4
	fstp	QWORD PTR [ebp-888]
L893:
	inc	edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L964
	fxch	st(1)
	fxch	st(2)
	fxch	st(3)
	fxch	st(4)
	xor	eax, eax
	cmp	eax, DWORD PTR [ebp+24]
	jl	L965
L950:
	fxch	st(4)
	fxch	st(3)
	fxch	st(2)
	fxch	st(1)
	jmp	L893
L964:
	fld	QWORD PTR [ebp-904]
	movsd	QWORD PTR [ebp-48], xmm3
	movsd	QWORD PTR [ebp-56], xmm1
	movsd	QWORD PTR [ebp-96], xmm0
	fstp	QWORD PTR [ebp-32]
	fld	QWORD PTR [ebp-896]
	fxch	st(5)
	fstp	QWORD PTR [ebp-104]
	fxch	st(2)
	fstp	QWORD PTR [ebp-64]
	fxch	st(1)
	fstp	QWORD PTR [ebp-72]
	fxch	st(2)
	fstp	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp-888]
	fxch	st(1)
	fstp	QWORD PTR [ebp-80]
	fxch	st(1)
	fstp	QWORD PTR [ebp-88]
	fstp	QWORD PTR [ebp-112]
	fld	QWORD PTR [ebp-880]
	fstp	QWORD PTR [ebp-120]
	jmp	L958
L932:
	fchs
	fmul	QWORD PTR [ebp+56]
	fmul	QWORD PTR [ebp+56]
	fstp	QWORD PTR [esp]
	call	_expm1
	fld	DWORD PTR LC179
	fdivrp	st(1), st
	fstp	QWORD PTR [ebp-736]
	jmp	L846
L934:
	fld	QWORD PTR [ebp+56]
	fmul	st, st(0)
	fdivr	QWORD PTR [ebp-760]
	fadd	QWORD PTR [ebp-736]
	fstp	QWORD PTR [ebp-736]
	jmp	L803
L947:
	fstp	st(0)
	fstp	st(0)
L901:
	inc	ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jl	L844
	jmp	L958
	.section .rdata,"dr"
	.align 8
LC185:
	.long	-1670671864
	.long	1072962135
	.align 4
LC186:
	.long	-1082130432
	.text
	.align 2
	.p2align 4,,15
.globl __Z5zinv2iPKdS0_S0_S0_PdS1_8SMOOTHERddb
	.def	__Z5zinv2iPKdS0_S0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z5zinv2iPKdS0_S0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 972
	mov	eax, DWORD PTR [ebp+36]
	fld	QWORD PTR [ebp+40]
	fxch	st(1)
	mov	edi, DWORD PTR [ebp+28]
	cmp	eax, 1
	movzx	edx, BYTE PTR [ebp+56]
	fstp	QWORD PTR [ebp-736]
	je	L994
	jle	L1177
	cmp	eax, 2
	je	L1089
	fstp	st(0)
	cmp	eax, 3
	je	L1178
L1202:
	add	esp, 972
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L994:
	fld	st(0)
	xor	ecx, ecx
	fmul	st, st(1)
	cmp	ecx, DWORD PTR [ebp+8]
	fdivr	QWORD PTR [ebp-736]
	fstp	QWORD PTR [ebp-768]
	jge	L1203
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	esi, [ebp-688]
	fld	QWORD PTR _I+8
	fxch	st(1)
	fstp	QWORD PTR [ebp-912]
	fstp	QWORD PTR [ebp-920]
	.p2align 4,,15
L1039:
	lea	ebx, [ecx+1]
	mov	edx, ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jge	L1147
	movsd	QWORD PTR [ebp-960], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-960]
	fld	QWORD PTR [ebp-912]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+16]
	fmul	st, st(2)
	fxch	st(1)
	fstp	QWORD PTR [ebp-800]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-776]
	fstp	QWORD PTR [ebp-808]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+24]
	fstp	QWORD PTR [ebp-816]
	fld	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [ebp-824]
	jmp	L1038
	.p2align 4,,7
L1186:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-200]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-192]
	jae	L1179
L1204:
	fld	QWORD PTR [ebp-128]
	fld	QWORD PTR [ebp-136]
	fld	QWORD PTR [ebp-768]
	fxch	st(2)
	fchs
	fxch	st(2)
	fmul	st, st(1)
	fld	st(2)
	fst	QWORD PTR [ebp-272]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fst	QWORD PTR [ebp-256]
	fmul	QWORD PTR [ebp-768]
	fxch	st(1)
	fsubrp	st(3), st
	fxch	st(1)
	fst	QWORD PTR [ebp-280]
	fst	QWORD PTR [ebp-264]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-296]
	fst	QWORD PTR [ebp-248]
	fstp	QWORD PTR [ebp-200]
	faddp	st(1), st
	fst	QWORD PTR [ebp-288]
	fst	QWORD PTR [ebp-240]
L1199:
	fstp	QWORD PTR [ebp-192]
	fld	QWORD PTR [ebp-200]
	mov	eax, DWORD PTR [ebp+20]
	fld	QWORD PTR [ebp-192]
	fld	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+24]
	fld	st(0)
	fld	QWORD PTR [eax+edx*8]
	fxch	st(1)
	fmul	st, st(4)
	mov	eax, DWORD PTR [ebp+32]
	fld	st(1)
	fmul	st, st(4)
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(2)
	fmul	st, st(5)
	fxch	st(1)
	fsubrp	st(3), st
	fxch	st(2)
	fsubr	QWORD PTR [edi+ecx*8]
	fxch	st(1)
	faddp	st(2), st
	fxch	st(1)
	fsubr	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+ecx*8]
	fld	QWORD PTR [ebp-824]
	fxch	st(1)
	fstp	QWORD PTR [eax+ecx*8]
	fmul	st, st(1)
	fld	QWORD PTR [ebp-816]
	fxch	st(2)
	fmul	QWORD PTR [ebp-816]
	fxch	st(2)
	fmul	st, st(3)
	fxch	st(3)
	fmul	QWORD PTR [ebp-824]
	fxch	st(3)
	fsubrp	st(1), st
	fadd	QWORD PTR [edi+edx*8]
	fxch	st(1)
	faddp	st(2), st
	fxch	st(1)
	fadd	QWORD PTR [eax+edx*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L1185
L1038:
	fld	QWORD PTR [ebp-800]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-808]
	fld	QWORD PTR [ebp-912]
	fxch	st(2)
	fsub	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+16]
	fld	QWORD PTR [ebp-920]
	fxch	st(2)
	fsub	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(3)
	lea	eax, [ebp-696]
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-920]
	fxch	st(2)
	fsubrp	st(3), st
	fxch	st(1)
	fadd	QWORD PTR [ebp-776]
	fxch	st(1)
	fadd	st, st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-168]
	fstp	QWORD PTR [ebp-152]
	fld	st(0)
	fld	st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-160]
	fxch	st(2)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fst	QWORD PTR [ebp-144]
	fst	QWORD PTR [ebp-176]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	QWORD PTR [ebp-696]
	fxch	st(2)
	cmovbe	eax, esi
	fst	QWORD PTR [ebp-184]
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fst	QWORD PTR [ebp-128]
	fxch	st(2)
	fstp	QWORD PTR [ebp-688]
	fld	QWORD PTR [eax]
	fucomi	st, st(3)
	fld	st(0)
	jp	L1186
	jne	L1186
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fxch	st(1)
	fst	QWORD PTR [ebp-200]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-192]
	jb	L1204
L1179:
	fld	QWORD PTR [ebp-136]
	fld	QWORD PTR [ebp-128]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fld	st(4)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-232]
	fxch	st(1)
	fst	QWORD PTR [ebp-224]
	fxch	st(1)
	fst	QWORD PTR [ebp-216]
	fxch	st(1)
	fst	QWORD PTR [ebp-208]
	fxch	st(1)
	fstp	QWORD PTR [ebp-200]
	jmp	L1199
L1089:
	fmul	st, st(0)
	test	dl, dl
	fdivr	QWORD PTR LC185
	fst	QWORD PTR [ebp-792]
	jne	L1180
	fstp	st(0)
L1090:
	xor	esi, esi
	.p2align 4,,15
L1175:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L1202
	mov	edx, DWORD PTR [ebp+8]
	lea	ecx, [esi+1]
	mov	ebx, ecx
	mov	DWORD PTR [ebp-860], ecx
	cmp	ecx, edx
	jge	L1154
	fldz
	jmp	L1137
	.p2align 4,,7
L1189:
	fstp	st(0)
	fdiv	st(2), st
	fdiv	st(1), st
	fxch	st(2)
	fmul	st, st(0)
	fxch	st(1)
	fmul	st, st(0)
	faddp	st(1), st
	fsqrt
	fmulp	st(1), st
L1110:
	fucomi	st, st(1)
	fld	st(1)
	jp	L1174
	je	L1191
L1174:
	fld	QWORD PTR [ebp+48]
	fxch	st(1)
	fst	QWORD PTR [ebp-584]
	fst	QWORD PTR [ebp-576]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L1205
	fstp	st(0)
	fld	QWORD PTR [ebp-520]
	fld	QWORD PTR [ebp-512]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fxch	st(4)
	fmul	st, st(3)
	fxch	st(4)
	faddp	st(3), st
	fxch	st(3)
	fsubrp	st(1), st
	fxch	st(1)
	fdiv	st, st(2)
	fxch	st(1)
	fdivrp	st(2), st
	fst	QWORD PTR [ebp-616]
	fxch	st(1)
	fst	QWORD PTR [ebp-608]
	fxch	st(1)
	fst	QWORD PTR [ebp-600]
	fxch	st(1)
	fst	QWORD PTR [ebp-592]
	fxch	st(1)
	fstp	QWORD PTR [ebp-584]
L1200:
	fstp	QWORD PTR [ebp-576]
	fld	QWORD PTR [ebp-584]
	mov	eax, DWORD PTR [ebp+20]
	mov	ecx, DWORD PTR [ebp+24]
	fld	QWORD PTR [ebp-576]
	mov	edx, DWORD PTR [ebp+32]
	fld	QWORD PTR [eax+ebx*8]
	mov	eax, DWORD PTR [ebp+20]
	fld	QWORD PTR [ecx+ebx*8]
	fld	st(1)
	mov	ecx, DWORD PTR [ebp+24]
	fmul	st, st(3)
	fld	st(1)
	fmul	st, st(5)
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(3)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fld	QWORD PTR [ecx+esi*8]
	fxch	st(2)
	fsubr	QWORD PTR [edx+esi*8]
	fxch	st(1)
	fsubrp	st(3), st
	fld	st(1)
	fxch	st(3)
	fsubr	QWORD PTR [edi+esi*8]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(2)
	fmul	st, st(5)
	fxch	st(1)
	fstp	QWORD PTR [edx+esi*8]
	mov	edx, DWORD PTR [ebp+32]
	fld	QWORD PTR [eax+esi*8]
	fxch	st(3)
	fstp	QWORD PTR [edi+esi*8]
	fld	st(2)
	fmulp	st(5), st
	fxch	st(2)
	fmulp	st(3), st
	fsubp	st(3), st
	faddp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR [edi+ebx*8]
	fxch	st(1)
	fadd	QWORD PTR [edx+ebx*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+ebx*8]
	fstp	QWORD PTR [edx+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jge	L1188
L1137:
	fld	QWORD PTR _I
	mov	ecx, DWORD PTR [ebp+16]
	lea	edx, [ebp-720]
	fld	QWORD PTR _I+8
	mov	eax, DWORD PTR [ebp+12]
	fld	st(1)
	fld	QWORD PTR [ecx+ebx*8]
	fxch	st(3)
	fmul	st, st(4)
	fld	QWORD PTR [eax+ebx*8]
	fld	QWORD PTR [eax+esi*8]
	fxch	st(5)
	lea	eax, [ebp-728]
	fsubr	QWORD PTR [ecx+esi*8]
	fxch	st(5)
	fsubrp	st(1), st
	fld	st(3)
	fmul	st, st(6)
	fxch	st(3)
	fmul	st, st(5)
	fxch	st(5)
	fmulp	st(4), st
	fxch	st(4)
	fsubrp	st(2), st
	faddp	st(2), st
	fadd	st(2), st
	fst	QWORD PTR [ebp-552]
	fld	st(1)
	fld	st(3)
	fxch	st(3)
	fst	QWORD PTR [ebp-544]
	fxch	st(3)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fstp	QWORD PTR [ebp-536]
	fxch	st(2)
	fst	QWORD PTR [ebp-528]
	fxch	st(1)
	fucomi	st, st(2)
	fxch	st(2)
	fstp	QWORD PTR [ebp-720]
	cmovbe	eax, edx
	fst	QWORD PTR [ebp-560]
	fxch	st(2)
	fst	QWORD PTR [ebp-568]
	fst	QWORD PTR [ebp-520]
	fxch	st(2)
	fst	QWORD PTR [ebp-512]
	fxch	st(1)
	fstp	QWORD PTR [ebp-728]
	fld	QWORD PTR [eax]
	fucomi	st, st(3)
	fld	st(0)
	jp	L1189
	jne	L1189
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	jmp	L1110
	.p2align 4,,7
L1185:
	fstp	st(0)
L1147:
	cmp	ebx, DWORD PTR [ebp+8]
	mov	ecx, ebx
	jl	L1039
L1203:
	fstp	st(0)
	add	esp, 972
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 4,,7
L1191:
	fstp	st(0)
	fstp	st(0)
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L1137
L1188:
	fstp	st(0)
L1154:
	mov	esi, DWORD PTR [ebp-860]
	jmp	L1175
L1205:
	fstp	st(1)
	fld	QWORD PTR [ebp-792]
	fchs
	fmul	st, st(1)
	fmulp	st(1), st
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-952]
	call	_expm1
	fld	QWORD PTR [ebp-520]
	fld	QWORD PTR [ebp-512]
	fxch	st(2)
	fchs
	fld	st(1)
	fld	st(3)
	fld	QWORD PTR [ebp-952]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(1)
	fmul	st, st(5)
	fld	st(3)
	fmul	st, st(5)
	fxch	st(2)
	faddp	st(1), st
	fxch	st(3)
	fmul	st, st(5)
	fxch	st(4)
	fmul	st, st(2)
	fxch	st(5)
	fmul	st, st(2)
	fxch	st(5)
	fsubrp	st(4), st
	faddp	st(4), st
	fld	QWORD PTR [ebp-736]
	fxch	st(3)
	fdiv	st, st(2)
	fxch	st(4)
	fdivrp	st(2), st
	fxch	st(3)
	fst	QWORD PTR [ebp-656]
	fld	st(0)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(1)
	fst	QWORD PTR [ebp-640]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fmul	QWORD PTR [ebp-736]
	fxch	st(2)
	fst	QWORD PTR [ebp-664]
	fxch	st(1)
	fsubrp	st(3), st
	fst	QWORD PTR [ebp-648]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-680]
	fst	QWORD PTR [ebp-632]
	fstp	QWORD PTR [ebp-584]
	faddp	st(1), st
	fst	QWORD PTR [ebp-672]
	fst	QWORD PTR [ebp-624]
	jmp	L1200
L1177:
	test	eax, eax
	jne	L1203
	fmul	st, st(0)
	test	dl, dl
	fstp	QWORD PTR [ebp-784]
	jne	L1182
L1041:
	xor	ecx, ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jge	L1202
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	esi, [ebp-704]
	fld	QWORD PTR _I+8
	fxch	st(1)
	fstp	QWORD PTR [ebp-896]
	fstp	QWORD PTR [ebp-904]
L1088:
	lea	ebx, [ecx+1]
	mov	edx, ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jge	L1151
	movsd	QWORD PTR [ebp-960], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-960]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+16]
	fstp	QWORD PTR [ebp-872]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-880]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+24]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	fstp	QWORD PTR [ebp-888]
	jmp	L1087
	.p2align 4,,7
L1194:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fld	QWORD PTR [ebp+48]
	fxch	st(3)
	fst	QWORD PTR [ebp-376]
	fst	QWORD PTR [ebp-368]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jae	L1183
L1206:
	fld	QWORD PTR [ebp-304]
	fxch	st(3)
	fmul	st, st(0)
	fld	QWORD PTR [ebp-312]
	fxch	st(4)
	fchs
	fxch	st(1)
	fadd	QWORD PTR [ebp-784]
	fld	st(4)
	fld	st(2)
	fxch	st(6)
	fst	QWORD PTR [ebp-472]
	fld	st(2)
	fmul	st, st(3)
	fxch	st(4)
	fst	QWORD PTR [ebp-464]
	fxch	st(2)
	fmul	st, st(3)
	fxch	st(7)
	fmul	st, st(5)
	fxch	st(1)
	fst	QWORD PTR [ebp-456]
	fmul	st, st(5)
	fxch	st(4)
	fadd	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-448]
	fmulp	st(3), st
	faddp	st(6), st
	fld	QWORD PTR [ebp-736]
	fxch	st(2)
	fsubrp	st(3), st
	fdiv	st(5), st
	fdivp	st(2), st
	fld	st(1)
	fxch	st(1)
	fmul	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-480]
	fxch	st(1)
	fmul	st, st(3)
	fxch	st(1)
	fst	QWORD PTR [ebp-432]
	fmul	QWORD PTR [ebp-736]
	fxch	st(2)
	fsubrp	st(1), st
	fxch	st(4)
	fst	QWORD PTR [ebp-488]
	fst	QWORD PTR [ebp-440]
	fmul	st, st(2)
	fxch	st(4)
	fst	QWORD PTR [ebp-504]
	fst	QWORD PTR [ebp-424]
	fstp	QWORD PTR [ebp-376]
	faddp	st(3), st
	fxch	st(2)
	fst	QWORD PTR [ebp-496]
	fst	QWORD PTR [ebp-416]
L1201:
	fstp	QWORD PTR [ebp-368]
	fld	QWORD PTR [ebp-376]
	mov	eax, DWORD PTR [ebp+20]
	fld	QWORD PTR [ebp-368]
	fld	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+24]
	fld	st(0)
	fld	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(3)
	mov	eax, DWORD PTR [ebp+32]
	fld	st(2)
	fmul	st, st(5)
	fxch	st(2)
	fmul	st, st(5)
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(2), st
	fxch	st(1)
	fsubr	QWORD PTR [eax+ecx*8]
	fxch	st(2)
	fsubrp	st(1), st
	fsubr	QWORD PTR [edi+ecx*8]
	fxch	st(1)
	fstp	QWORD PTR [eax+ecx*8]
	fld	QWORD PTR [ebp-888]
	fxch	st(1)
	fstp	QWORD PTR [edi+ecx*8]
	fld	st(3)
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-888]
	fxch	st(1)
	fmul	st, st(3)
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fsubrp	st(2), st
	faddp	st(2), st
	fadd	QWORD PTR [edi+edx*8]
	fxch	st(1)
	fadd	QWORD PTR [eax+edx*8]
	fxch	st(1)
	fstp	QWORD PTR [edi+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L1193
L1087:
	fld	QWORD PTR [ebp-872]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-880]
	fld	QWORD PTR [ebp-896]
	fxch	st(2)
	fsub	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+16]
	fld	QWORD PTR [ebp-904]
	fld	QWORD PTR [ebp-896]
	fxch	st(3)
	fsub	QWORD PTR [eax+edx*8]
	fxch	st(1)
	fmul	st, st(6)
	fxch	st(3)
	lea	eax, [ebp-712]
	fmul	st, st(6)
	fxch	st(4)
	fmul	st, st(1)
	fxch	st(1)
	fmul	QWORD PTR [ebp-904]
	fxch	st(1)
	fsubrp	st(3), st
	faddp	st(3), st
	fxch	st(1)
	fst	QWORD PTR [ebp-344]
	fadd	st(1), st
	fld	st(2)
	fst	QWORD PTR [ebp-336]
	fld	st(2)
	fabs
	fxch	st(2)
	fstp	QWORD PTR [ebp-328]
	fxch	st(3)
	fabs
	fxch	st(3)
	fst	QWORD PTR [ebp-320]
	fxch	st(1)
	fst	QWORD PTR [ebp-704]
	fxch	st(3)
	fucomi	st, st(3)
	fstp	st(3)
	cmovbe	eax, esi
	fst	QWORD PTR [ebp-352]
	fxch	st(1)
	fst	QWORD PTR [ebp-360]
	fst	QWORD PTR [ebp-312]
	fxch	st(1)
	fst	QWORD PTR [ebp-304]
	fxch	st(2)
	fstp	QWORD PTR [ebp-712]
	fld	QWORD PTR [eax]
	fucomi	st, st(4)
	fld	st(0)
	jp	L1194
	jne	L1194
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fld	QWORD PTR [ebp+48]
	fxch	st(3)
	fst	QWORD PTR [ebp-376]
	fst	QWORD PTR [ebp-368]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jb	L1206
L1183:
	fstp	st(2)
	fld	QWORD PTR [ebp-312]
	fld	QWORD PTR [ebp-304]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(6)
	fxch	st(1)
	faddp	st(2), st
	fld	st(5)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-408]
	fxch	st(1)
	fst	QWORD PTR [ebp-400]
	fxch	st(1)
	fst	QWORD PTR [ebp-392]
	fxch	st(1)
	fst	QWORD PTR [ebp-384]
	fxch	st(1)
	fstp	QWORD PTR [ebp-376]
	jmp	L1201
L1178:
	xor	edx, edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L1202
	fld	QWORD PTR _I
	fld	QWORD PTR _I+8
	fld	QWORD PTR [ebp-88]
	fld	QWORD PTR [ebp-32]
	fxch	st(3)
	fstp	QWORD PTR [ebp-928]
	fxch	st(1)
	fstp	QWORD PTR [ebp-936]
	fld	QWORD PTR [ebp-104]
	fld	QWORD PTR [ebp-80]
	fld	QWORD PTR [ebp-40]
	movsd	xmm3, QWORD PTR [ebp-112]
	movsd	xmm2, QWORD PTR [ebp-120]
	movsd	xmm1, QWORD PTR [ebp-48]
	movsd	xmm0, QWORD PTR [ebp-56]
	fld	QWORD PTR [ebp-96]
	fld	QWORD PTR [ebp-64]
	fld	QWORD PTR [ebp-72]
	lea	ecx, [edx+1]
	mov	eax, ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jge	L1196
L1208:
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	mov	esi, DWORD PTR [ebp+12]
	fldz
	mov	ebx, DWORD PTR [ebp+16]
	fmul	QWORD PTR [ebp-936]
	fld	QWORD PTR [esi+edx*8]
	mov	esi, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-832]
	fld	QWORD PTR [ebx+edx*8]
	fxch	st(1)
	mov	ebx, DWORD PTR [ebp+24]
	fstp	QWORD PTR [ebp-752]
	fstp	QWORD PTR [ebp-840]
	fld	QWORD PTR [esi+edx*8]
	fstp	QWORD PTR [ebp-848]
	fld	QWORD PTR [ebx+edx*8]
	fstp	QWORD PTR [ebp-856]
	fldz
	fmul	QWORD PTR [ebp-928]
	fstp	QWORD PTR [ebp-760]
	jmp	L992
	.p2align 4,,7
L1197:
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
	fstp	st(0)
L992:
	fld	QWORD PTR [ebp-840]
	mov	esi, DWORD PTR [ebp+16]
	fldz
	fld	QWORD PTR [ebp-928]
	mov	ebx, DWORD PTR [ebp+12]
	fldz
	fxch	st(3)
	fsub	QWORD PTR [esi+eax*8]
	mov	esi, DWORD PTR [ebp+24]
	fld	QWORD PTR [ebp-936]
	fld	QWORD PTR [ebp-832]
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(1)
	fmulp	st(2), st
	fxch	st(2)
	fsub	QWORD PTR [ebx+eax*8]
	fxch	st(2)
	mov	ebx, DWORD PTR [ebp+20]
	fsub	QWORD PTR [ebp-752]
	fxch	st(1)
	fadd	QWORD PTR [ebp-760]
	fld	QWORD PTR [esi+eax*8]
	fxch	st(3)
	fadd	st, st(2)
	fld	st(1)
	fxch	st(3)
	fstp	QWORD PTR [ebp-744]
	fld	st(0)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(3)
	fmul	st, st(1)
	fxch	st(5)
	fmul	st, st(2)
	fxch	st(6)
	fmul	st, st(1)
	fxch	st(5)
	faddp	st(3), st
	fadd	st(5), st
	fld	st(3)
	fld	QWORD PTR [ebx+eax*8]
	fxch	st(6)
	fsub	st, st(3)
	fxch	st(7)
	mov	ebx, DWORD PTR [ebp+32]
	fdiv	st, st(4)
	fxch	st(7)
	fdivrp	st(4), st
	fld	st(5)
	fmul	st, st(7)
	fxch	st(1)
	fmul	st, st(4)
	fxch	st(6)
	fmul	st, st(4)
	fxch	st(5)
	fmul	st, st(7)
	fxch	st(1)
	fsubrp	st(6), st
	fxch	st(5)
	fsubr	QWORD PTR [edi+edx*8]
	fxch	st(4)
	faddp	st(5), st
	fxch	st(4)
	fsubr	QWORD PTR [ebx+edx*8]
	fxch	st(3)
	fstp	QWORD PTR [edi+edx*8]
	fld	QWORD PTR [ebp-856]
	fxch	st(3)
	fstp	QWORD PTR [ebx+edx*8]
	fxch	st(2)
	fmul	st, st(1)
	fld	QWORD PTR [ebp-848]
	fmul	st, st(5)
	fsubrp	st(1), st
	fld	QWORD PTR [ebp-856]
	fxch	st(1)
	fadd	QWORD PTR [edi+eax*8]
	fxch	st(1)
	fmul	st, st(5)
	fxch	st(1)
	fstp	QWORD PTR [edi+eax*8]
	fld	QWORD PTR [ebp-848]
	fmul	st, st(2)
	faddp	st(1), st
	fadd	QWORD PTR [ebx+eax*8]
	fstp	QWORD PTR [ebx+eax*8]
	inc	eax
	cmp	eax, DWORD PTR [ebp+8]
	jl	L1197
	fst	QWORD PTR [ebp-960]
	movsd	xmm3, QWORD PTR [ebp-960]
	fld	st(2)
	fxch	st(4)
	fst	QWORD PTR [ebp-960]
	fld	st(2)
	movsd	xmm2, QWORD PTR [ebp-960]
	fst	QWORD PTR [ebp-960]
	fld	QWORD PTR [ebp-744]
	fxch	st(5)
	movsd	xmm1, QWORD PTR [ebp-960]
	fstp	QWORD PTR [ebp-960]
	fld	st(4)
	movsd	xmm0, QWORD PTR [ebp-960]
	fld	st(1)
L1143:
	cmp	ecx, DWORD PTR [ebp+8]
	mov	edx, ecx
	jge	L1207
	fxch	st(2)
	fxch	st(4)
	fxch	st(7)
	fxch	st(5)
	fxch	st(1)
	lea	ecx, [edx+1]
	mov	eax, ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jl	L1208
L1196:
	fxch	st(1)
	fxch	st(5)
	fxch	st(7)
	fxch	st(4)
	fxch	st(2)
	jmp	L1143
L1207:
	fxch	st(1)
	fstp	QWORD PTR [ebp-72]
	fxch	st(4)
	fstp	QWORD PTR [ebp-64]
	fxch	st(3)
	fstp	QWORD PTR [ebp-96]
	movsd	QWORD PTR [ebp-56], xmm0
	movsd	QWORD PTR [ebp-48], xmm1
	movsd	QWORD PTR [ebp-120], xmm2
	movsd	QWORD PTR [ebp-112], xmm3
	fstp	QWORD PTR [ebp-40]
	fstp	QWORD PTR [ebp-32]
	fxch	st(2)
	fstp	QWORD PTR [ebp-104]
	fxch	st(1)
	fstp	QWORD PTR [ebp-80]
	fstp	QWORD PTR [ebp-88]
	jmp	L1202
L1180:
	fchs
	fmul	QWORD PTR [ebp+48]
	fmul	QWORD PTR [ebp+48]
	fstp	QWORD PTR [esp]
	call	_expm1
	fld	DWORD PTR LC186
	fdivrp	st(1), st
	fstp	QWORD PTR [ebp-736]
	jmp	L1090
L1182:
	fld	QWORD PTR [ebp+48]
	fmul	st, st(0)
	fdivr	QWORD PTR [ebp-784]
	fadd	QWORD PTR [ebp-736]
	fstp	QWORD PTR [ebp-736]
	jmp	L1041
L1193:
	fstp	st(0)
	fstp	st(0)
L1151:
	cmp	ebx, DWORD PTR [ebp+8]
	mov	ecx, ebx
	jl	L1088
	jmp	L1202
	.section .rdata,"dr"
	.align 8
LC192:
	.long	-1670671864
	.long	1072962135
	.align 4
LC193:
	.long	-1082130432
	.text
	.align 2
	.p2align 4,,15
.globl __Z5rinv2iPKdS0_S0_PdS1_8SMOOTHERddb
	.def	__Z5rinv2iPKdS0_S0_PdS1_8SMOOTHERddb;	.scl	2;	.type	32;	.endef
__Z5rinv2iPKdS0_S0_PdS1_8SMOOTHERddb:
	push	ebp
	fld1
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 860
	fld	QWORD PTR [ebp+36]
	fxch	st(1)
	movzx	ecx, BYTE PTR [ebp+52]
	fstp	QWORD PTR [ebp-648]
	mov	BYTE PTR [ebp-633], cl
	mov	DWORD PTR [esp], 176
	fstp	QWORD PTR [ebp-840]
	call	_mxMalloc
	cmp	DWORD PTR [ebp+32], 1
	mov	ecx, eax
	lea	edx, [eax+16]
	mov	DWORD PTR [ebp-652], edx
	lea	eax, [eax+32]
	lea	edi, [ecx+48]
	mov	DWORD PTR [ebp-656], eax
	lea	esi, [ecx+64]
	lea	ebx, [ecx+80]
	mov	DWORD PTR [ebp-660], edi
	lea	edx, [ecx+144]
	lea	eax, [ecx+160]
	mov	DWORD PTR [ebp-664], esi
	lea	edi, [ecx+96]
	lea	esi, [ecx+112]
	mov	DWORD PTR [ebp-668], ebx
	lea	ebx, [ecx+128]
	mov	DWORD PTR [ebp-672], edx
	mov	DWORD PTR [ebp-676], eax
	fld	QWORD PTR [ebp-840]
	je	L1238
	jle	L1407
	cmp	DWORD PTR [ebp+32], 2
	je	L1325
	fstp	st(0)
	cmp	DWORD PTR [ebp+32], 3
	je	L1408
L1431:
	add	esp, 860
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
L1238:
	fld	st(0)
	xor	ecx, ecx
	fmul	st, st(1)
	cmp	ecx, DWORD PTR [ebp+8]
	fdivr	QWORD PTR [ebp-648]
	fstp	QWORD PTR [ebp-688]
	jge	L1432
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	esi, [ebp-592]
	fld	QWORD PTR _I+8
	fxch	st(1)
	lea	edi, [ebp-600]
	fstp	QWORD PTR [ebp-792]
	fstp	QWORD PTR [ebp-800]
	.p2align 4,,15
L1279:
	lea	ebx, [ecx+1]
	mov	edx, ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jge	L1379
	movsd	QWORD PTR [ebp-848], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-848]
	fld	QWORD PTR [ebp-792]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+16]
	fmul	st, st(2)
	fxch	st(1)
	fstp	QWORD PTR [ebp-728]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-696]
	fstp	QWORD PTR [ebp-736]
	fld	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [ebp-744]
	jmp	L1278
	.p2align 4,,7
L1417:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-104]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-96]
	jae	L1409
L1433:
	fld	QWORD PTR [ebp-32]
	fld	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp-688]
	fxch	st(2)
	fchs
	fxch	st(2)
	fmul	st, st(1)
	fld	st(2)
	fst	QWORD PTR [ebp-176]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fst	QWORD PTR [ebp-160]
	fmul	QWORD PTR [ebp-688]
	fxch	st(1)
	fsubrp	st(3), st
	fxch	st(1)
	fst	QWORD PTR [ebp-184]
	fst	QWORD PTR [ebp-168]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-200]
	fst	QWORD PTR [ebp-152]
	fstp	QWORD PTR [ebp-104]
	faddp	st(1), st
	fst	QWORD PTR [ebp-192]
	fst	QWORD PTR [ebp-144]
L1427:
	fstp	QWORD PTR [ebp-96]
	fld	QWORD PTR [ebp-104]
	mov	eax, DWORD PTR [ebp+20]
	fld	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+24]
	fld	st(0)
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-744]
	fxch	st(2)
	fsubr	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp-96]
	fmul	st(1), st
	fmul	QWORD PTR [ebp-744]
	fxch	st(1)
	fsubr	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+24]
	fadd	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+28]
	fadd	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L1416
L1278:
	fld	QWORD PTR [ebp-728]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-736]
	fld	QWORD PTR [ebp-792]
	fxch	st(2)
	fsub	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+16]
	fld	QWORD PTR [ebp-800]
	fxch	st(2)
	fsub	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(3)
	mov	eax, edi
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-800]
	fxch	st(2)
	fsubrp	st(3), st
	fxch	st(1)
	fadd	QWORD PTR [ebp-696]
	fxch	st(1)
	fadd	st, st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-72]
	fstp	QWORD PTR [ebp-56]
	fld	st(0)
	fld	st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-64]
	fxch	st(2)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fst	QWORD PTR [ebp-48]
	fst	QWORD PTR [ebp-80]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	QWORD PTR [ebp-600]
	fxch	st(2)
	cmovbe	eax, esi
	fst	QWORD PTR [ebp-88]
	fst	QWORD PTR [ebp-40]
	fxch	st(1)
	fst	QWORD PTR [ebp-32]
	fxch	st(2)
	fstp	QWORD PTR [ebp-592]
	fld	QWORD PTR [eax]
	fucomi	st, st(3)
	fld	st(0)
	jp	L1417
	jne	L1417
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fxch	st(1)
	fst	QWORD PTR [ebp-104]
	fxch	st(1)
	fucomip	st, st(2)
	fst	QWORD PTR [ebp-96]
	jb	L1433
L1409:
	fld	QWORD PTR [ebp-40]
	fld	QWORD PTR [ebp-32]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fld	st(4)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-136]
	fxch	st(1)
	fst	QWORD PTR [ebp-128]
	fxch	st(1)
	fst	QWORD PTR [ebp-120]
	fxch	st(1)
	fst	QWORD PTR [ebp-112]
	fxch	st(1)
	fstp	QWORD PTR [ebp-104]
	jmp	L1427
L1325:
	cmp	BYTE PTR [ebp-633], 0
	fmul	st, st(0)
	fdivr	QWORD PTR LC192
	fst	QWORD PTR [ebp-720]
	jne	L1410
	fstp	st(0)
L1326:
	xor	esi, esi
	.p2align 4,,15
L1405:
	cmp	esi, DWORD PTR [ebp+8]
	jge	L1431
	lea	edi, [esi+1]
	mov	ebx, edi
	cmp	edi, DWORD PTR [ebp+8]
	jge	L1386
	fldz
	jmp	L1369
	.p2align 4,,7
L1420:
	fstp	st(0)
	fdiv	st(2), st
	fdiv	st(1), st
	fxch	st(2)
	fmul	st, st(0)
	fxch	st(1)
	fmul	st, st(0)
	faddp	st(1), st
	fsqrt
	fmulp	st(1), st
L1346:
	fucomi	st, st(1)
	fld	st(1)
	jp	L1404
	je	L1422
L1404:
	fld	QWORD PTR [ebp+44]
	fxch	st(1)
	fst	QWORD PTR [ebp-488]
	fst	QWORD PTR [ebp-480]
	fxch	st(2)
	fucomi	st, st(1)
	fstp	st(1)
	jb	L1434
	fstp	st(0)
	fld	QWORD PTR [ebp-424]
	fld	QWORD PTR [ebp-416]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(5)
	fxch	st(1)
	faddp	st(2), st
	fxch	st(4)
	fmul	st, st(3)
	fxch	st(4)
	faddp	st(3), st
	fxch	st(3)
	fsubrp	st(1), st
	fxch	st(1)
	fdiv	st, st(2)
	fxch	st(1)
	fdivrp	st(2), st
	fst	QWORD PTR [ebp-520]
	fxch	st(1)
	fst	QWORD PTR [ebp-512]
	fxch	st(1)
	fst	QWORD PTR [ebp-504]
	fxch	st(1)
	fst	QWORD PTR [ebp-496]
	fxch	st(1)
	fstp	QWORD PTR [ebp-488]
L1428:
	fstp	QWORD PTR [ebp-480]
	fld	QWORD PTR [ebp-488]
	mov	ecx, DWORD PTR [ebp+20]
	mov	eax, DWORD PTR [ebp+24]
	mov	edx, DWORD PTR [ebp+28]
	fld	QWORD PTR [ecx+ebx*8]
	mov	ecx, DWORD PTR [ebp+20]
	fld	st(0)
	fmul	st, st(2)
	fsubr	QWORD PTR [eax+esi*8]
	fstp	QWORD PTR [eax+esi*8]
	mov	eax, DWORD PTR [ebp+24]
	fld	QWORD PTR [ebp-480]
	fld	QWORD PTR [ecx+esi*8]
	fxch	st(2)
	fmul	st, st(1)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(3)
	fsubr	QWORD PTR [edx+esi*8]
	fxch	st(2)
	fmulp	st(1), st
	fxch	st(2)
	fadd	QWORD PTR [eax+ebx*8]
	fxch	st(1)
	fstp	QWORD PTR [edx+esi*8]
	mov	edx, DWORD PTR [ebp+28]
	fstp	QWORD PTR [eax+ebx*8]
	fadd	QWORD PTR [edx+ebx*8]
	fstp	QWORD PTR [edx+ebx*8]
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jge	L1419
L1369:
	fld	QWORD PTR _I
	mov	eax, DWORD PTR [ebp+16]
	lea	edx, [ebp-624]
	fld	QWORD PTR _I+8
	mov	ecx, DWORD PTR [ebp+12]
	fld	st(1)
	fld	QWORD PTR [eax+ebx*8]
	fxch	st(3)
	fmul	st, st(4)
	fld	QWORD PTR [ecx+ebx*8]
	fld	QWORD PTR [ecx+esi*8]
	fxch	st(5)
	lea	ecx, [ebp-632]
	fsubr	QWORD PTR [eax+esi*8]
	fxch	st(5)
	fsubrp	st(1), st
	fld	st(3)
	fmul	st, st(6)
	fxch	st(3)
	fmul	st, st(5)
	fxch	st(5)
	fmulp	st(4), st
	fxch	st(4)
	fsubrp	st(2), st
	faddp	st(2), st
	fadd	st(2), st
	fst	QWORD PTR [ebp-456]
	fld	st(1)
	fld	st(3)
	fxch	st(3)
	fst	QWORD PTR [ebp-448]
	fxch	st(3)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fstp	QWORD PTR [ebp-440]
	fxch	st(2)
	fst	QWORD PTR [ebp-432]
	fxch	st(1)
	fucomi	st, st(2)
	fxch	st(2)
	fstp	QWORD PTR [ebp-624]
	cmovbe	ecx, edx
	fst	QWORD PTR [ebp-464]
	fxch	st(2)
	fst	QWORD PTR [ebp-472]
	fst	QWORD PTR [ebp-424]
	fxch	st(2)
	fst	QWORD PTR [ebp-416]
	fxch	st(1)
	fstp	QWORD PTR [ebp-632]
	fld	QWORD PTR [ecx]
	fucomi	st, st(3)
	fld	st(0)
	jp	L1420
	jne	L1420
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	jmp	L1346
L1408:
	xor	eax, eax
	mov	edx, DWORD PTR [ebp+8]
	mov	DWORD PTR [ebp-680], eax
	cmp	DWORD PTR [ebp-680], edx
	jge	L1431
	fld	QWORD PTR [ebp-648]
	lea	eax, [edx-1]
	mov	DWORD PTR [ebp-764], eax
L1237:
	mov	DWORD PTR [ecx], 0
	mov	edx, DWORD PTR [ebp-680]
	mov	eax, DWORD PTR [ebp+12]
	mov	DWORD PTR [ecx+4], 1072693248
	mov	DWORD PTR [ecx+8], 0
	mov	DWORD PTR [ecx+12], 1072693248
	fld	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [ecx+16]
	fld	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+16]
	fstp	QWORD PTR [ecx+24]
	fld	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [ecx+32]
	fld	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+20]
	mov	DWORD PTR [ecx+48], 0
	mov	DWORD PTR [ecx+52], 0
	fstp	QWORD PTR [ecx+40]
	mov	DWORD PTR [ecx+56], 0
	mov	DWORD PTR [ecx+60], 0
	mov	DWORD PTR [ecx+64], 0
	mov	DWORD PTR [ecx+68], 0
	mov	DWORD PTR [ecx+72], 0
	mov	DWORD PTR [ecx+76], 0
	fld	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [ecx+80]
	fld	QWORD PTR [eax+edx*8]
	inc	edx
	test	dl, 1
	mov	DWORD PTR [ebp-768], edx
	fstp	QWORD PTR [ecx+88]
	mov	DWORD PTR [ebp-852], edx
	je	L1215
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [eax+edx*8]
	mov	edx, DWORD PTR [ebp-680]
	fsubr	QWORD PTR [eax+edx*8]
	mov	edx, DWORD PTR [ebp+16]
	mov	eax, DWORD PTR [ebp-768]
	fld	st(0)
	fmul	st, st(1)
	fld	QWORD PTR [edx+eax*8]
	mov	eax, DWORD PTR [ebp-680]
	fsubr	QWORD PTR [edx+eax*8]
	mov	edx, DWORD PTR [ebp-768]
	mov	eax, DWORD PTR [ebp+20]
	fld	st(0)
	fmul	st, st(1)
	fld	QWORD PTR [eax+edx*8]
	fxch	st(3)
	mov	edx, DWORD PTR [ebp-680]
	mov	eax, DWORD PTR [ebp+24]
	faddp	st(1), st
	fdivr	st, st(4)
	fmul	st(2), st
	fld	st(2)
	fmul	st, st(4)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(3)
	fsubr	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	fxch	st(2)
	mov	eax, DWORD PTR [ebp+28]
	fadd	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+20]
	fmul	QWORD PTR [eax+edx*8]
	mov	edx, DWORD PTR [ebp-768]
	mov	eax, DWORD PTR [ebp+24]
	fmul	st(2), st
	fmulp	st(1), st
	fxch	st(1)
	fadd	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+28]
	fsubr	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	mov	edx, DWORD PTR [ebp-680]
	add	edx, 2
	mov	DWORD PTR [ebp-852], edx
	jmp	L1215
	.p2align 4,,7
L1414:
	mov	eax, DWORD PTR [ebp-652]
	mov	edx, DWORD PTR [ebp+12]
	movapd	xmm6, XMMWORD PTR [eax]
	mov	eax, DWORD PTR [ebp-852]
	sar	eax
	sal	eax, 4
	addpd	xmm6, XMMWORD PTR [eax+edx]
	mov	edx, DWORD PTR [ebp-656]
	movapd	XMMWORD PTR [edi], xmm6
	movapd	xmm5, XMMWORD PTR [edx]
	mov	edx, DWORD PTR [ebp+16]
	addpd	xmm5, XMMWORD PTR [eax+edx]
	mov	edx, DWORD PTR [ebp-672]
	movapd	XMMWORD PTR [esi], xmm5
	movapd	xmm4, XMMWORD PTR [edi]
	mulpd	xmm4, xmm4
	movapd	XMMWORD PTR [ebx], xmm4
	movapd	xmm3, XMMWORD PTR [esi]
	mulpd	xmm3, xmm3
	movapd	XMMWORD PTR [edx], xmm3
	movapd	xmm1, XMMWORD PTR [ebx]
	addpd	xmm1, xmm3
	movapd	XMMWORD PTR [edx], xmm1
	mov	edx, DWORD PTR [ebp-676]
	movapd	xmm2, XMMWORD PTR [ecx]
	divpd	xmm2, xmm1
	movapd	XMMWORD PTR [edx], xmm2
	mulpd	xmm2, XMMWORD PTR [edi]
	movapd	XMMWORD PTR [edi], xmm2
	movapd	xmm0, XMMWORD PTR [edx]
	mov	edx, DWORD PTR [ebp+20]
	mulpd	xmm0, XMMWORD PTR [esi]
	movapd	XMMWORD PTR [esi], xmm0
	movapd	xmm7, XMMWORD PTR [eax+edx]
	mulpd	xmm7, XMMWORD PTR [edi]
	movapd	XMMWORD PTR [ebx], xmm7
	movapd	xmm6, XMMWORD PTR [eax+edx]
	mov	edx, DWORD PTR [ebp-672]
	mulpd	xmm6, XMMWORD PTR [esi]
	movapd	XMMWORD PTR [edx], xmm6
	mov	edx, DWORD PTR [ebp-660]
	add	DWORD PTR [ebp-852], 2
	movapd	xmm5, XMMWORD PTR [edx]
	subpd	xmm5, XMMWORD PTR [ebx]
	movapd	XMMWORD PTR [edx], xmm5
	mov	edx, DWORD PTR [ebp-664]
	movapd	xmm4, XMMWORD PTR [edx]
	mov	edx, DWORD PTR [ebp-672]
	addpd	xmm4, XMMWORD PTR [edx]
	mov	edx, DWORD PTR [ebp-664]
	movapd	XMMWORD PTR [edx], xmm4
	mov	edx, DWORD PTR [ebp-668]
	movapd	xmm3, XMMWORD PTR [edx]
	mulpd	xmm3, XMMWORD PTR [edi]
	movapd	XMMWORD PTR [ebx], xmm3
	movapd	xmm2, XMMWORD PTR [edx]
	mov	edx, DWORD PTR [ebp-672]
	mulpd	xmm2, XMMWORD PTR [esi]
	movapd	XMMWORD PTR [edx], xmm2
	mov	edx, DWORD PTR [ebp+24]
	movapd	xmm1, XMMWORD PTR [eax+edx]
	addpd	xmm1, XMMWORD PTR [ebx]
	movapd	XMMWORD PTR [eax+edx], xmm1
	mov	edx, DWORD PTR [ebp+28]
	movapd	xmm0, XMMWORD PTR [eax+edx]
	subpd	xmm0, XMMWORD PTR [ebx]
	movapd	XMMWORD PTR [eax+edx], xmm0
L1215:
	mov	eax, DWORD PTR [ebp-852]
	cmp	DWORD PTR [ebp-764], eax
	jg	L1414
	fld	QWORD PTR [ecx+56]
	mov	edx, DWORD PTR [ebp-680]
	mov	eax, DWORD PTR [ebp+24]
	fadd	QWORD PTR [ecx+48]
	fadd	QWORD PTR [eax+edx*8]
	fst	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+28]
	fld	QWORD PTR [ecx+72]
	fadd	QWORD PTR [ecx+64]
	fadd	QWORD PTR [eax+edx*8]
	mov	edx, DWORD PTR [ebp+8]
	cmp	DWORD PTR [ebp-852], edx
	jge	L1387
	mov	edx, DWORD PTR [ebp+12]
	mov	eax, DWORD PTR [ebp-852]
	fld	QWORD PTR [edx+eax*8]
	mov	eax, DWORD PTR [ebp-680]
	fsubr	QWORD PTR [edx+eax*8]
	mov	eax, DWORD PTR [ebp+16]
	mov	edx, DWORD PTR [ebp-852]
	fld	st(0)
	fmul	st, st(1)
	fld	QWORD PTR [eax+edx*8]
	mov	edx, DWORD PTR [ebp-680]
	fsubr	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp-852]
	mov	edx, DWORD PTR [ebp+20]
	fld	st(0)
	fmul	st, st(1)
	faddp	st(2), st
	fld	QWORD PTR [edx+eax*8]
	fxch	st(2)
	mov	eax, DWORD PTR [ebp-680]
	fdivr	st, st(6)
	mov	edx, DWORD PTR [ebp+24]
	fmul	st(2), st
	fld	st(2)
	fmul	st, st(4)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(6)
	fsubrp	st(3), st
	fxch	st(4)
	faddp	st(5), st
	fxch	st(1)
	fstp	QWORD PTR [edx+eax*8]
	fxch	st(3)
	mov	edx, DWORD PTR [ebp+28]
	fstp	QWORD PTR [edx+eax*8]
	fxch	st(1)
	mov	edx, DWORD PTR [ebp+20]
	fmul	QWORD PTR [edx+eax*8]
	mov	eax, DWORD PTR [ebp-852]
	mov	edx, DWORD PTR [ebp+24]
	fmul	st(1), st
	fmulp	st(2), st
	fadd	QWORD PTR [edx+eax*8]
	fstp	QWORD PTR [edx+eax*8]
	mov	edx, DWORD PTR [ebp+28]
	fsubr	QWORD PTR [edx+eax*8]
L1430:
	fstp	QWORD PTR [edx+eax*8]
	mov	eax, DWORD PTR [ebp-768]
	mov	edx, DWORD PTR [ebp+8]
	mov	DWORD PTR [ebp-680], eax
	cmp	eax, edx
	jl	L1237
L1432:
	fstp	st(0)
L1435:
	add	esp, 860
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 4,,7
L1422:
	fstp	st(0)
	fstp	st(0)
	inc	ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jl	L1369
L1419:
	fstp	st(0)
L1386:
	mov	esi, edi
	jmp	L1405
	.p2align 4,,7
L1416:
	fstp	st(0)
L1379:
	cmp	ebx, DWORD PTR [ebp+8]
	mov	ecx, ebx
	jl	L1279
	fstp	st(0)
	jmp	L1435
L1434:
	fstp	st(1)
	fld	QWORD PTR [ebp-720]
	fchs
	fmul	st, st(1)
	fmulp	st(1), st
	fstp	QWORD PTR [esp]
	fstp	QWORD PTR [ebp-824]
	call	_expm1
	fld	QWORD PTR [ebp-424]
	fld	QWORD PTR [ebp-416]
	fxch	st(2)
	fchs
	fld	st(1)
	fld	st(3)
	fld	QWORD PTR [ebp-824]
	fxch	st(2)
	fmul	st, st(4)
	fxch	st(1)
	fmul	st, st(5)
	fld	st(3)
	fmul	st, st(5)
	fxch	st(2)
	faddp	st(1), st
	fxch	st(3)
	fmul	st, st(5)
	fxch	st(4)
	fmul	st, st(2)
	fxch	st(5)
	fmul	st, st(2)
	fxch	st(5)
	fsubrp	st(4), st
	faddp	st(4), st
	fld	QWORD PTR [ebp-648]
	fxch	st(3)
	fdiv	st, st(2)
	fxch	st(4)
	fdivrp	st(2), st
	fxch	st(3)
	fst	QWORD PTR [ebp-560]
	fld	st(0)
	fxch	st(3)
	fmul	st, st(2)
	fxch	st(1)
	fst	QWORD PTR [ebp-544]
	fxch	st(3)
	fmul	st, st(4)
	fxch	st(3)
	fmul	QWORD PTR [ebp-648]
	fxch	st(2)
	fst	QWORD PTR [ebp-568]
	fxch	st(1)
	fsubrp	st(3), st
	fst	QWORD PTR [ebp-552]
	fmul	st, st(3)
	fxch	st(2)
	fst	QWORD PTR [ebp-584]
	fst	QWORD PTR [ebp-536]
	fstp	QWORD PTR [ebp-488]
	faddp	st(1), st
	fst	QWORD PTR [ebp-576]
	fst	QWORD PTR [ebp-528]
	jmp	L1428
L1407:
	mov	ebx, DWORD PTR [ebp+32]
	test	ebx, ebx
	jne	L1432
	cmp	BYTE PTR [ebp-633], 0
	fmul	st, st(0)
	fstp	QWORD PTR [ebp-704]
	jne	L1412
L1281:
	xor	ecx, ecx
	cmp	ecx, DWORD PTR [ebp+8]
	jge	L1431
	fld	QWORD PTR _I
	pxor	xmm0, xmm0
	lea	esi, [ebp-608]
	fld	QWORD PTR _I+8
	fxch	st(1)
	lea	edi, [ebp-616]
	fstp	QWORD PTR [ebp-776]
	fstp	QWORD PTR [ebp-784]
L1324:
	lea	ebx, [ecx+1]
	mov	edx, ebx
	cmp	ebx, DWORD PTR [ebp+8]
	jge	L1383
	movsd	QWORD PTR [ebp-848], xmm0
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-848]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+16]
	fstp	QWORD PTR [ebp-752]
	fld	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+20]
	fstp	QWORD PTR [ebp-760]
	fld	QWORD PTR [ebp-776]
	fld	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	fmul	st, st(2)
	fstp	QWORD PTR [ebp-712]
	jmp	L1323
	.p2align 4,,7
L1425:
	fstp	st(0)
	fdiv	st(1), st
	fdiv	st(2), st
	fxch	st(1)
	fmul	st, st(0)
	fxch	st(2)
	fmul	st, st(0)
	faddp	st(2), st
	fxch	st(1)
	fsqrt
	fmulp	st(1), st
	fld	QWORD PTR [ebp+44]
	fxch	st(3)
	fst	QWORD PTR [ebp-280]
	fst	QWORD PTR [ebp-272]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jae	L1413
L1436:
	fld	QWORD PTR [ebp-208]
	fxch	st(3)
	fmul	st, st(0)
	fld	QWORD PTR [ebp-216]
	fxch	st(4)
	fchs
	fxch	st(1)
	fadd	QWORD PTR [ebp-704]
	fld	st(4)
	fld	st(2)
	fxch	st(6)
	fst	QWORD PTR [ebp-376]
	fld	st(2)
	fmul	st, st(3)
	fxch	st(4)
	fst	QWORD PTR [ebp-368]
	fxch	st(2)
	fmul	st, st(3)
	fxch	st(7)
	fmul	st, st(5)
	fxch	st(1)
	fst	QWORD PTR [ebp-360]
	fmul	st, st(5)
	fxch	st(4)
	fadd	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-352]
	fmulp	st(3), st
	faddp	st(6), st
	fld	QWORD PTR [ebp-648]
	fxch	st(2)
	fsubrp	st(3), st
	fdiv	st(5), st
	fdivp	st(2), st
	fld	st(1)
	fxch	st(1)
	fmul	st, st(5)
	fxch	st(2)
	fst	QWORD PTR [ebp-384]
	fxch	st(1)
	fmul	st, st(3)
	fxch	st(1)
	fst	QWORD PTR [ebp-336]
	fmul	QWORD PTR [ebp-648]
	fxch	st(2)
	fsubrp	st(1), st
	fxch	st(4)
	fst	QWORD PTR [ebp-392]
	fst	QWORD PTR [ebp-344]
	fmul	st, st(2)
	fxch	st(4)
	fst	QWORD PTR [ebp-408]
	fst	QWORD PTR [ebp-328]
	fstp	QWORD PTR [ebp-280]
	faddp	st(3), st
	fxch	st(2)
	fst	QWORD PTR [ebp-400]
	fst	QWORD PTR [ebp-320]
L1429:
	fstp	QWORD PTR [ebp-272]
	fld	QWORD PTR [ebp-280]
	mov	eax, DWORD PTR [ebp+20]
	fld	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+24]
	fld	st(0)
	fmul	st, st(2)
	fxch	st(2)
	fmul	st, st(3)
	fxch	st(2)
	fsubr	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [eax+ecx*8]
	mov	eax, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp-272]
	fmul	st(1), st
	fmul	st, st(3)
	fxch	st(1)
	fsubr	QWORD PTR [eax+ecx*8]
	fstp	QWORD PTR [eax+ecx*8]
	fxch	st(1)
	mov	eax, DWORD PTR [ebp+24]
	fadd	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+28]
	fadd	QWORD PTR [eax+edx*8]
	fstp	QWORD PTR [eax+edx*8]
	inc	edx
	cmp	edx, DWORD PTR [ebp+8]
	jge	L1424
L1323:
	fld	QWORD PTR [ebp-752]
	mov	eax, DWORD PTR [ebp+12]
	fld	QWORD PTR [ebp-760]
	fld	QWORD PTR [ebp-776]
	fxch	st(2)
	fsub	QWORD PTR [eax+edx*8]
	mov	eax, DWORD PTR [ebp+16]
	fld	QWORD PTR [ebp-784]
	fxch	st(2)
	fsub	QWORD PTR [eax+edx*8]
	fxch	st(2)
	fmul	st, st(5)
	fxch	st(3)
	mov	eax, edi
	fmul	st, st(2)
	fxch	st(2)
	fmul	QWORD PTR [ebp-784]
	fxch	st(2)
	fsubrp	st(3), st
	fxch	st(1)
	fadd	QWORD PTR [ebp-712]
	fxch	st(1)
	fadd	st, st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-248]
	fstp	QWORD PTR [ebp-232]
	fld	st(0)
	fld	st(2)
	fxch	st(2)
	fst	QWORD PTR [ebp-240]
	fxch	st(2)
	fabs
	fxch	st(1)
	fabs
	fxch	st(2)
	fst	QWORD PTR [ebp-224]
	fst	QWORD PTR [ebp-256]
	fxch	st(2)
	fucomi	st, st(1)
	fxch	st(1)
	fstp	QWORD PTR [ebp-608]
	fxch	st(2)
	cmovbe	eax, esi
	fst	QWORD PTR [ebp-264]
	fst	QWORD PTR [ebp-216]
	fxch	st(1)
	fst	QWORD PTR [ebp-208]
	fxch	st(2)
	fstp	QWORD PTR [ebp-616]
	fld	QWORD PTR [eax]
	fucomi	st, st(4)
	fld	st(0)
	jp	L1425
	jne	L1425
	fstp	st(1)
	fstp	st(1)
	fstp	st(1)
	fld	QWORD PTR [ebp+44]
	fxch	st(3)
	fst	QWORD PTR [ebp-280]
	fst	QWORD PTR [ebp-272]
	fxch	st(1)
	fucomi	st, st(3)
	fstp	st(3)
	jb	L1436
L1413:
	fstp	st(2)
	fld	QWORD PTR [ebp-216]
	fld	QWORD PTR [ebp-208]
	fld	st(1)
	fld	st(1)
	fmul	st, st(2)
	fxch	st(1)
	fmul	st, st(3)
	fld	st(2)
	fmul	st, st(6)
	fxch	st(1)
	faddp	st(2), st
	fld	st(5)
	fmul	st, st(4)
	fxch	st(1)
	faddp	st(4), st
	fsubrp	st(2), st
	fdiv	st(2), st
	fdivp	st(1), st
	fxch	st(1)
	fst	QWORD PTR [ebp-312]
	fxch	st(1)
	fst	QWORD PTR [ebp-304]
	fxch	st(1)
	fst	QWORD PTR [ebp-296]
	fxch	st(1)
	fst	QWORD PTR [ebp-288]
	fxch	st(1)
	fstp	QWORD PTR [ebp-280]
	jmp	L1429
L1387:
	fstp	st(1)
	mov	eax, DWORD PTR [ebp-680]
	mov	edx, DWORD PTR [ebp+28]
	jmp	L1430
L1410:
	fchs
	fmul	QWORD PTR [ebp+44]
	fmul	QWORD PTR [ebp+44]
	fstp	QWORD PTR [esp]
	call	_expm1
	fld	DWORD PTR LC193
	fdivrp	st(1), st
	fstp	QWORD PTR [ebp-648]
	jmp	L1326
L1412:
	fld	QWORD PTR [ebp+44]
	fmul	st, st(0)
	fdivr	QWORD PTR [ebp-704]
	fadd	QWORD PTR [ebp-648]
	fstp	QWORD PTR [ebp-648]
	jmp	L1281
L1424:
	fstp	st(0)
	fstp	st(0)
L1383:
	cmp	ebx, DWORD PTR [ebp+8]
	mov	ecx, ebx
	jl	L1324
	jmp	L1431
	.section .rdata,"dr"
	.align 4
LC198:
	.ascii "Direct summation System time: %f\12\0"
	.align 4
LC197:
	.long	1148846080
	.text
	.align 2
	.p2align 4,,15
.globl __Z14directInteractiPKdS0_S0_S0_iS0_S0_PdS1_S1_S1_PK5panelii8SMOOTHERddbS1_i
	.def	__Z14directInteractiPKdS0_S0_S0_iS0_S0_PdS1_S1_S1_PK5panelii8SMOOTHERddbS1_i;	.scl	2;	.type	32;	.endef
__Z14directInteractiPKdS0_S0_S0_iS0_S0_PdS1_S1_S1_PK5panelii8SMOOTHERddbS1_i:
	push	ebp
	mov	ebp, esp
	push	edi
	push	esi
	push	ebx
	sub	esp, 76
	mov	eax, DWORD PTR [ebp+96]
	movzx	ebx, BYTE PTR [ebp+88]
	mov	ecx, DWORD PTR [ebp+92]
	test	eax, eax
	setne	dl
	mov	esi, DWORD PTR [ebp+40]
	test	ecx, ecx
	setne	al
	mov	edi, DWORD PTR [ebp+68]
	or	dl, al
	jne	L1456
L1438:
	test	esi, esi
	je	L1439
	mov	edx, DWORD PTR [ebp+24]
	test	edx, edx
	je	L1457
	mov	eax, DWORD PTR [ebp+64]
	test	eax, eax
	jne	L1458
	fld	QWORD PTR [ebp+80]
	movzx	ecx, bl
	mov	DWORD PTR [esp+28], edi
	mov	edx, DWORD PTR [ebp+44]
	mov	DWORD PTR [esp+48], ecx
	mov	eax, DWORD PTR [ebp+24]
	mov	DWORD PTR [esp+20], esi
	mov	ecx, DWORD PTR [ebp+20]
	fstp	QWORD PTR [esp+40]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+24], edx
	mov	edx, DWORD PTR [ebp+16]
	mov	DWORD PTR [esp+16], eax
	mov	eax, DWORD PTR [ebp+12]
	mov	DWORD PTR [esp+12], ecx
	mov	ecx, DWORD PTR [ebp+8]
	fstp	QWORD PTR [esp+32]
	mov	DWORD PTR [esp+8], edx
	mov	DWORD PTR [esp+4], eax
	mov	DWORD PTR [esp], ecx
	call	__Z5zlog2iPKdS0_S0_S0_PdS1_8SMOOTHERddb
L1443:
	mov	DWORD PTR [esp+4], esi
	mov	edx, DWORD PTR [ebp+8]
	mov	eax, DWORD PTR [ebp+60]
	mov	ecx, DWORD PTR [ebp+16]
	mov	esi, DWORD PTR [ebp+56]
	mov	DWORD PTR [esp+24], edx
	mov	edx, DWORD PTR [ebp+12]
	mov	DWORD PTR [esp+20], eax
	mov	eax, DWORD PTR [ebp+44]
	mov	DWORD PTR [esp+16], ecx
	mov	DWORD PTR [esp+12], edx
	mov	DWORD PTR [esp+8], eax
	mov	DWORD PTR [esp], esi
	call	__Z19directInteractPanelPK5panelPdS2_PKdS4_ii
L1439:
	mov	esi, DWORD PTR [ebp+48]
	test	esi, esi
	je	L1446
	mov	ecx, DWORD PTR [ebp+24]
	test	ecx, ecx
	je	L1459
	mov	esi, DWORD PTR [ebp+64]
	test	esi, esi
	je	L1451
	fld	QWORD PTR [ebp+80]
	movzx	edx, bl
	mov	DWORD PTR [esp+40], edi
	mov	ebx, DWORD PTR [ebp+52]
	mov	DWORD PTR [esp+60], edx
	mov	eax, DWORD PTR [ebp+48]
	mov	ecx, DWORD PTR [ebp+36]
	fstp	QWORD PTR [esp+52]
	mov	edi, DWORD PTR [ebp+32]
	mov	esi, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+36], ebx
	mov	edx, DWORD PTR [ebp+24]
	mov	DWORD PTR [esp+32], eax
	mov	ebx, DWORD PTR [ebp+20]
	mov	DWORD PTR [esp+28], ecx
	mov	eax, DWORD PTR [ebp+16]
	mov	DWORD PTR [esp+24], edi
	mov	ecx, DWORD PTR [ebp+12]
	fstp	QWORD PTR [esp+44]
	mov	edi, DWORD PTR [ebp+8]
	mov	DWORD PTR [esp+20], esi
	mov	DWORD PTR [esp+16], edx
	mov	DWORD PTR [esp+12], ebx
	mov	DWORD PTR [esp+8], eax
	mov	DWORD PTR [esp+4], ecx
	mov	DWORD PTR [esp], edi
	call	__Z4zinviPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb
L1450:
	mov	ebx, DWORD PTR [ebp+28]
	mov	edx, DWORD PTR [ebp+60]
	mov	eax, DWORD PTR [ebp+36]
	mov	DWORD PTR [esp+24], ebx
	mov	ecx, DWORD PTR [ebp+32]
	mov	edi, DWORD PTR [ebp+52]
	mov	DWORD PTR [esp+20], edx
	mov	esi, DWORD PTR [ebp+48]
	mov	ebx, DWORD PTR [ebp+56]
	mov	DWORD PTR [esp+16], eax
	mov	DWORD PTR [esp+12], ecx
	mov	DWORD PTR [esp+8], edi
	mov	DWORD PTR [esp+4], esi
	mov	DWORD PTR [esp], ebx
	call	__Z19directInteractPanelPK5panelPdS2_PKdS4_ii
L1446:
	mov	edi, DWORD PTR [ebp+96]
	mov	esi, DWORD PTR [ebp+92]
	test	edi, edi
	setne	dl
	test	esi, esi
	setne	cl
	or	dl, cl
	jne	L1460
	mov	eax, DWORD PTR [ebp+96]
	test	eax, eax
	jne	L1461
L1454:
	mov	edx, DWORD PTR [ebp+92]
	test	edx, edx
	je	L1437
	mov	edi, DWORD PTR _TIME_after
	mov	eax, DWORD PTR _TIME_before
	mov	esi, DWORD PTR [ebp+92]
	sub	edi, eax
	cvtsi2sd	xmm1, edi
	movsd	QWORD PTR [ebp-24], xmm1
	fld	QWORD PTR [ebp-24]
	fdiv	DWORD PTR LC197
	fstp	QWORD PTR [esi]
L1437:
	add	esp, 76
	pop	ebx
	pop	esi
	pop	edi
	pop	ebp
	ret
	.p2align 4,,7
L1458:
	fld	QWORD PTR [ebp+80]
	movzx	edx, bl
	mov	DWORD PTR [esp+28], edi
	mov	eax, DWORD PTR [ebp+44]
	mov	DWORD PTR [esp+48], edx
	mov	ecx, DWORD PTR [ebp+24]
	mov	DWORD PTR [esp+20], esi
	mov	edx, DWORD PTR [ebp+20]
	fstp	QWORD PTR [esp+40]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+24], eax
	mov	eax, DWORD PTR [ebp+16]
	mov	DWORD PTR [esp+16], ecx
	mov	ecx, DWORD PTR [ebp+12]
	mov	DWORD PTR [esp+12], edx
	mov	edx, DWORD PTR [ebp+8]
	fstp	QWORD PTR [esp+32]
	mov	DWORD PTR [esp+8], eax
	mov	DWORD PTR [esp+4], ecx
	mov	DWORD PTR [esp], edx
	call	__Z5zinv2iPKdS0_S0_S0_PdS1_8SMOOTHERddb
	jmp	L1443
	.p2align 4,,7
L1451:
	fld	QWORD PTR [ebp+80]
	mov	DWORD PTR [esp+40], edi
	movzx	edx, bl
	mov	eax, DWORD PTR [ebp+52]
	mov	DWORD PTR [esp+60], edx
	mov	ecx, DWORD PTR [ebp+48]
	mov	edi, DWORD PTR [ebp+36]
	fstp	QWORD PTR [esp+52]
	mov	esi, DWORD PTR [ebp+32]
	mov	edx, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+36], eax
	mov	ebx, DWORD PTR [ebp+24]
	mov	DWORD PTR [esp+32], ecx
	mov	eax, DWORD PTR [ebp+20]
	mov	DWORD PTR [esp+28], edi
	mov	ecx, DWORD PTR [ebp+16]
	mov	DWORD PTR [esp+24], esi
	mov	edi, DWORD PTR [ebp+12]
	fstp	QWORD PTR [esp+44]
	mov	esi, DWORD PTR [ebp+8]
	mov	DWORD PTR [esp+20], edx
	mov	DWORD PTR [esp+16], ebx
	mov	DWORD PTR [esp+12], eax
	mov	DWORD PTR [esp+8], ecx
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], esi
	call	__Z4zlogiPKdS0_S0_S0_iS0_S0_PdS1_8SMOOTHERddb
	jmp	L1450
	.p2align 4,,7
L1456:
	call	_clock
	mov	DWORD PTR _TIME_before, eax
	jmp	L1438
	.p2align 4,,7
L1460:
	call	_clock
	mov	DWORD PTR _TIME_after, eax
	mov	eax, DWORD PTR [ebp+96]
	test	eax, eax
	je	L1454
	.p2align 4,,15
L1461:
	mov	DWORD PTR [esp], OFFSET FLAT:LC198
	mov	ecx, DWORD PTR _TIME_before
	mov	ebx, DWORD PTR _TIME_after
	sub	ebx, ecx
	cvtsi2sd	xmm0, ebx
	movsd	QWORD PTR [ebp-24], xmm0
	fld	QWORD PTR [ebp-24]
	fdiv	DWORD PTR LC197
	fstp	QWORD PTR [esp+4]
	call	_mexPrintf
	jmp	L1454
	.p2align 4,,7
L1459:
	mov	eax, DWORD PTR [ebp+64]
	test	eax, eax
	je	L1448
	fld	QWORD PTR [ebp+80]
	movzx	esi, bl
	mov	DWORD PTR [esp+36], edi
	mov	edx, DWORD PTR [ebp+52]
	mov	DWORD PTR [esp+56], esi
	mov	edi, DWORD PTR [ebp+48]
	mov	ebx, DWORD PTR [ebp+36]
	fstp	QWORD PTR [esp+48]
	mov	eax, DWORD PTR [ebp+32]
	mov	ecx, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+32], edx
	mov	esi, DWORD PTR [ebp+20]
	mov	DWORD PTR [esp+28], edi
	mov	edx, DWORD PTR [ebp+16]
	mov	DWORD PTR [esp+24], ebx
	mov	edi, DWORD PTR [ebp+12]
	mov	DWORD PTR [esp+20], eax
	mov	ebx, DWORD PTR [ebp+8]
	fstp	QWORD PTR [esp+40]
	mov	DWORD PTR [esp+16], ecx
	mov	DWORD PTR [esp+12], esi
	mov	DWORD PTR [esp+8], edx
	mov	DWORD PTR [esp+4], edi
	mov	DWORD PTR [esp], ebx
	call	__Z4rinviPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb
	jmp	L1450
	.p2align 4,,7
L1457:
	mov	ecx, DWORD PTR [ebp+64]
	test	ecx, ecx
	je	L1441
	fld	QWORD PTR [ebp+80]
	movzx	ecx, bl
	mov	DWORD PTR [esp+24], edi
	mov	edx, DWORD PTR [ebp+44]
	mov	DWORD PTR [esp+44], ecx
	mov	eax, DWORD PTR [ebp+20]
	mov	DWORD PTR [esp+16], esi
	mov	ecx, DWORD PTR [ebp+16]
	fstp	QWORD PTR [esp+36]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+20], edx
	mov	edx, DWORD PTR [ebp+12]
	mov	DWORD PTR [esp+12], eax
	mov	eax, DWORD PTR [ebp+8]
	fstp	QWORD PTR [esp+28]
	mov	DWORD PTR [esp+8], ecx
	mov	DWORD PTR [esp+4], edx
	mov	DWORD PTR [esp], eax
	call	__Z5rinv2iPKdS0_S0_PdS1_8SMOOTHERddb
	jmp	L1443
L1448:
	fld	QWORD PTR [ebp+80]
	movzx	edx, bl
	mov	DWORD PTR [esp+36], edi
	mov	ebx, DWORD PTR [ebp+52]
	mov	DWORD PTR [esp+56], edx
	mov	eax, DWORD PTR [ebp+48]
	mov	ecx, DWORD PTR [ebp+36]
	fstp	QWORD PTR [esp+48]
	mov	esi, DWORD PTR [ebp+32]
	mov	edx, DWORD PTR [ebp+28]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+32], ebx
	mov	edi, DWORD PTR [ebp+20]
	mov	DWORD PTR [esp+28], eax
	mov	ebx, DWORD PTR [ebp+16]
	mov	DWORD PTR [esp+24], ecx
	mov	eax, DWORD PTR [ebp+12]
	mov	DWORD PTR [esp+20], esi
	mov	ecx, DWORD PTR [ebp+8]
	fstp	QWORD PTR [esp+40]
	mov	DWORD PTR [esp+16], edx
	mov	DWORD PTR [esp+12], edi
	mov	DWORD PTR [esp+8], ebx
	mov	DWORD PTR [esp+4], eax
	mov	DWORD PTR [esp], ecx
	call	__Z4rlogiPKdS0_S0_iS0_S0_PdS1_8SMOOTHERddb
	jmp	L1450
L1441:
	fld	QWORD PTR [ebp+80]
	movzx	ecx, bl
	mov	DWORD PTR [esp+24], edi
	mov	edx, DWORD PTR [ebp+44]
	mov	DWORD PTR [esp+44], ecx
	mov	eax, DWORD PTR [ebp+20]
	mov	DWORD PTR [esp+16], esi
	mov	ecx, DWORD PTR [ebp+16]
	fstp	QWORD PTR [esp+36]
	fld	QWORD PTR [ebp+72]
	mov	DWORD PTR [esp+20], edx
	mov	edx, DWORD PTR [ebp+12]
	mov	DWORD PTR [esp+12], eax
	mov	eax, DWORD PTR [ebp+8]
	fstp	QWORD PTR [esp+28]
	mov	DWORD PTR [esp+8], ecx
	mov	DWORD PTR [esp+4], edx
	mov	DWORD PTR [esp], eax
	call	__Z5rlog2iPKdS0_S0_PdS1_8SMOOTHERddb
	jmp	L1443
	.align 2
	.p2align 4,,15
	.def	__GLOBAL__I__Z14directInteractiPKdS0_S0_S0_iS0_S0_PdS1_S1_S1_PK5panelii8SMOOTHERddbS1_i;	.scl	3;	.type	32;	.endef
__GLOBAL__I__Z14directInteractiPKdS0_S0_S0_iS0_S0_PdS1_S1_S1_PK5panelii8SMOOTHERddbS1_i:
	push	ebp
	fldz
	mov	ebp, esp
	pop	ebp
	fstp	QWORD PTR _I
	fld1
	fstp	QWORD PTR _I+8
	ret
	.def	_expm1;	.scl	2;	.type	32;	.endef
	.def	_mxMalloc;	.scl	2;	.type	32;	.endef
	.def	_mexPrintf;	.scl	2;	.type	32;	.endef
	.def	__Z19directInteractPanelPK5panelPdS2_PKdS4_ii;	.scl	2;	.type	32;	.endef
	.def	_clock;	.scl	2;	.type	32;	.endef
	.def	_log;	.scl	2;	.type	32;	.endef
	.def	_hypot;	.scl	2;	.type	32;	.endef
	.def	_exp;	.scl	2;	.type	32;	.endef
