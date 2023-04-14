main:
 stp	x29, x30, [sp, #-64]!
 adrp	x1, 411000 <__libc_start_main@GLIBC_2.17>
 adrp	x0, 9ca6000 <A+0x9894e88>
 add	x3, x1, #0x178
 mov	x29, sp
 fmov	s2, #4.000000000000000000e+00
 fmov	s5, #-4.000000000000000000e+00
 str	x23, [sp, #48]
 stp	x19, x20, [sp, #16]
 stp	x21, x22, [sp, #32]
 stp	s5, s2, [x3, #4]
 add	x0, x0, #0xd78
 add	x2, x0, #0xc28
 adrp	x21, 411000 <__libc_start_main@GLIBC_2.17>
 fmov	s6, #3.000000000000000000e+00
 adrp	x20, 400000 <_init-0x6d0>
 fmov	s4, #5.000000000000000000e+00
 add	x21, x21, #0x60
 fmov	s3, #1.000000000000000000e+00
 add	x20, x20, #0xb10
 fmov	s0, #-3.000000000000000000e+00
 fmov	s1, #2.000000000000000000e+00
 str	s6, [x1, #376]
 adrp	x1, 1353d000 <B+0x9895660>
 str	s2, [x0, #3112]
 add	x1, x1, #0xa78
 str	s3, [x0, #3116]
 add	x19, x1, #0x750
 str	s4, [x0, #3120]
 add	x23, x1, #0x760
 str	s0, [x3, #12]
 ldr	q0, [x3]
 str	s1, [x0, #3124]
 ldr	q5, [x2]
 fmov	s3, #6.000000000000000000e+00
 adrp	x0, 411000 <__libc_start_main@GLIBC_2.17>
 str	s2, [x2, #28]
 add	x0, x0, #0x50
 stp	s3, s6, [x2, #16]
 fmov	s6, #-2.000000000000000000e+00
 trn2	v1.4s, v0.4s, v0.4s
 ldr	q7, [x0]
 trn1	v0.4s, v0.4s, v0.4s
 fmov	s2, #7.000000000000000000e+00
 rev64	v16.4s, v5.4s
 str	s3, [x3, #24]
 stp	s4, s6, [x3, #16]
 fmov	s4, #-1.000000000000000000e+00
 str	s2, [x2, #24]
 fmul	v0.4s, v0.4s, v5.4s
 str	s4, [x3, #28]
 fmul	v1.4s, v1.4s, v16.4s
 fmla	v0.4s, v1.4s, v7.4s
 str	q0, [x1, #1872]
 ldr	s0, [x19]
 mov	x0, x21
 add	x19, x19, #0x8
 fcvt	d0, s0
 bl	400780 <std::ostream& std::ostream::_M_insert<double>(double)@plt>
 mov	x22, x0
 mov	x2, #0x1                   	// #1
 mov	x1, x20
 bl	400770 <std::basic_ostream<char, std::char_traits<char> >& std::__ostream_insert<char, std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*, long)@plt>
 ldur	s0, [x19, #-4]
 mov	x0, x22
 fcvt	d0, s0
 bl	400780 <std::ostream& std::ostream::_M_insert<double>(double)@plt>
 mov	x2, #0x1                   	// #1
 mov	x1, x20
 bl	400770 <std::basic_ostream<char, std::char_traits<char> >& std::__ostream_insert<char, std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*, long)@plt>
 cmp	x19, x23
 b.ne	400868 <main+0xd8>  // b.any
 mov	w0, #0x0                   	// #0
 ldr	x23, [sp, #48]
 ldp	x19, x20, [sp, #16]
 ldp	x21, x22, [sp, #32]
 ldp	x29, x30, [sp], #64
 ret
_GLOBAL__sub_I_A:
 stp	x29, x30, [sp, #-32]!
 mov	x29, sp
 str	x19, [sp, #16]
 adrp	x19, 1cdd4000 <C+0x9895e38>
 add	x19, x19, #0x778
 add	x19, x19, #0x278
 mov	x0, x19
 bl	400730 <std::ios_base::Init::Init()@plt>
 mov	x1, x19
 adrp	x2, 411000 <__libc_start_main@GLIBC_2.17>
 ldr	x19, [sp, #16]
 adrp	x0, 400000 <_init-0x6d0>
 ldp	x29, x30, [sp], #32
 add	x2, x2, #0x48
 add	x0, x0, #0x750
 b	400740 <__cxa_atexit@plt>
call_weak_fn:
 adrp	x0, 410000 <__FRAME_END__+0xf430>
 ldr	x0, [x0, #4064]
 cbz	x0, 400960 <call_weak_fn+0x10>
 b	400720 <__gmon_start__@plt>
 ret
 udf	#0