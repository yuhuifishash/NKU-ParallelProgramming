call_weak_fn:
 adrp	x0, 411000 <__FRAME_END__+0xff30>
 ldr	x0, [x0, #4064]
 cbz	x0, 400798 <call_weak_fn+0x10>
 b	4006d0 <__gmon_start__@plt>
 ret
 udf	#0
V_complex_mul(Complex*, Complex*, Complex*, float*):
 sub	sp, sp, #0x240
 str	x0, [sp, #72]
 str	x1, [sp, #64]
 str	x2, [sp, #56]
 str	x3, [sp, #48]
 ldr	x0, [sp, #72]
 str	x0, [sp, #360]
 ldr	x0, [sp, #360]
 ldr	q0, [x0]
 str	q0, [sp, #544]
 ldr	x0, [sp, #64]
 str	x0, [sp, #568]
 ldr	x0, [sp, #568]
 ldr	q0, [x0]
 str	q0, [sp, #512]
 ldr	q0, [sp, #544]
 str	q0, [sp, #464]
 ldr	q0, [sp, #544]
 str	q0, [sp, #256]
 ldr	q0, [sp, #464]
 str	q0, [sp, #240]
 ldr	q0, [sp, #256]
 str	q0, [sp, #224]
 ldr	q0, 400d10 <V_complex_mul(Complex*, Complex*, Complex*, float*)+0x490>
 mov	w0, v0.s[0]
 and	w0, w0, #0x4
 mov	w1, v0.s[0]
 and	w2, w1, #0x3
 ldr	q3, [sp, #240]
 str	q3, [sp, #112]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x3, sp, #0x240
 add	x1, x3, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s4, [x1, #3632]
 ldr	q3, [sp, #224]
 str	q3, [sp, #128]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3648]
 cmp	w0, wzr
 fcsel	s3, s4, s3, eq  // eq = none
 fmov	w5, s3
 mov	w0, v0.s[1]
 and	w0, w0, #0x4
 mov	w1, v0.s[1]
 and	w2, w1, #0x3
 mov	w1, w2
 lsl	x1, x1, #2
 add	x3, sp, #0x240
 add	x1, x3, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s4, [x1, #3632]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3648]
 cmp	w0, wzr
 fcsel	s3, s4, s3, eq  // eq = none
 fmov	w4, s3
 mov	w0, v0.s[2]
 and	w0, w0, #0x4
 mov	w1, v0.s[2]
 and	w2, w1, #0x3
 mov	w1, w2
 lsl	x1, x1, #2
 add	x3, sp, #0x240
 add	x1, x3, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s4, [x1, #3632]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3648]
 cmp	w0, wzr
 fcsel	s3, s4, s3, eq  // eq = none
 fmov	w3, s3
 mov	w0, v0.s[3]
 and	w0, w0, #0x4
 mov	w1, v0.s[3]
 and	w2, w1, #0x3
 mov	w1, w2
 lsl	x1, x1, #2
 add	x6, sp, #0x240
 add	x1, x6, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3632]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s0, [x1, #3648]
 cmp	w0, wzr
 fcsel	s0, s3, s0, eq  // eq = none
 fmov	w0, s0
 str	w5, [sp, #32]
 str	w4, [sp, #36]
 str	w3, [sp, #40]
 str	w0, [sp, #44]
 ldr	q0, [sp, #32]
 mov	v1.16b, v0.16b
 ldr	q0, [sp, #464]
 str	q0, [sp, #208]
 ldr	q0, [sp, #256]
 str	q0, [sp, #192]
 ldr	q0, 400d20 <V_complex_mul(Complex*, Complex*, Complex*, float*)+0x4a0>
 mov	w0, v0.s[0]
 and	w0, w0, #0x4
 mov	w1, v0.s[0]
 and	w2, w1, #0x3
 ldr	q3, [sp, #208]
 str	q3, [sp, #144]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x3, sp, #0x240
 add	x1, x3, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s4, [x1, #3664]
 ldr	q3, [sp, #192]
 str	q3, [sp, #160]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3680]
 cmp	w0, wzr
 fcsel	s3, s4, s3, eq  // eq = none
 fmov	w5, s3
 mov	w0, v0.s[1]
 and	w0, w0, #0x4
 mov	w1, v0.s[1]
 and	w2, w1, #0x3
 mov	w1, w2
 lsl	x1, x1, #2
 add	x3, sp, #0x240
 add	x1, x3, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s4, [x1, #3664]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3680]
 cmp	w0, wzr
 fcsel	s3, s4, s3, eq  // eq = none
 fmov	w4, s3
 mov	w0, v0.s[2]
 and	w0, w0, #0x4
 mov	w1, v0.s[2]
 and	w2, w1, #0x3
 mov	w1, w2
 lsl	x1, x1, #2
 add	x3, sp, #0x240
 add	x1, x3, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s4, [x1, #3664]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3680]
 cmp	w0, wzr
 fcsel	s3, s4, s3, eq  // eq = none
 fmov	w3, s3
 mov	w0, v0.s[3]
 and	w0, w0, #0x4
 mov	w1, v0.s[3]
 and	w2, w1, #0x3
 mov	w1, w2
 lsl	x1, x1, #2
 add	x6, sp, #0x240
 add	x1, x6, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s3, [x1, #3664]
 mov	w1, w2
 lsl	x1, x1, #2
 add	x2, sp, #0x240
 add	x1, x2, x1
 sub	x1, x1, #0x1, lsl #12
 ldr	s0, [x1, #3680]
 cmp	w0, wzr
 fcsel	s0, s3, s0, eq  // eq = none
 fmov	w0, s0
 str	w5, [sp, #16]
 str	w4, [sp, #20]
 str	w3, [sp, #24]
 str	w0, [sp, #28]
 ldr	q0, [sp, #16]
 mov	v2.16b, v0.16b
 add	x0, sp, #0x50
 st1	{v1.16b, v2.16b}, [x0]
 add	x0, sp, #0x50
 ldr	q0, [x0, #16]
 str	q0, [sp, #496]
 add	x0, sp, #0x50
 ldr	q0, [x0]
 str	q0, [sp, #480]
 ldr	q0, [sp, #512]
 str	q0, [sp, #432]
 ldr	q1, [sp, #432]
 ldr	q0, 400d30 <V_complex_mul(Complex*, Complex*, Complex*, float*)+0x4b0>
 mov	w0, v0.s[0]
 and	w0, w0, #0x3
 str	q1, [sp, #176]
 mov	w0, w0
 lsl	x0, x0, #2
 add	x1, sp, #0x240
 add	x0, x1, x0
 sub	x0, x0, #0x1, lsl #12
 ldr	w3, [x0, #3696]
 mov	w0, v0.s[1]
 and	w0, w0, #0x3
 mov	w0, w0
 lsl	x0, x0, #2
 add	x1, sp, #0x240
 add	x0, x1, x0
 sub	x0, x0, #0x1, lsl #12
 ldr	w2, [x0, #3696]
 mov	w0, v0.s[2]
 and	w0, w0, #0x3
 mov	w0, w0
 lsl	x0, x0, #2
 add	x1, sp, #0x240
 add	x0, x1, x0
 sub	x0, x0, #0x1, lsl #12
 ldr	w1, [x0, #3696]
 mov	w0, v0.s[3]
 and	w0, w0, #0x3
 mov	w0, w0
 lsl	x0, x0, #2
 add	x4, sp, #0x240
 add	x0, x4, x0
 sub	x0, x0, #0x1, lsl #12
 ldr	w0, [x0, #3696]
 str	w3, [sp]
 str	w2, [sp, #4]
 str	w1, [sp, #8]
 str	w0, [sp, #12]
 ldr	q0, [sp]
 str	q0, [sp, #448]
 ldr	q0, [sp, #496]
 str	q0, [sp, #400]
 ldr	q0, [sp, #448]
 str	q0, [sp, #272]
 ldr	q1, [sp, #400]
 ldr	q0, [sp, #272]
 fmul	v0.4s, v1.4s, v0.4s
 str	q0, [sp, #416]
 ldr	x0, [sp, #48]
 ldr	q0, [x0]
 ldr	q1, [sp, #416]
 str	q1, [sp, #384]
 str	q0, [sp, #288]
 ldr	q1, [sp, #384]
 ldr	q0, [sp, #288]
 fmul	v0.4s, v1.4s, v0.4s
 str	q0, [sp, #416]
 ldr	q0, [sp, #416]
 str	q0, [sp, #368]
 ldr	q0, [sp, #480]
 str	q0, [sp, #320]
 ldr	q0, [sp, #512]
 str	q0, [sp, #304]
 ldr	q1, [sp, #320]
 ldr	q0, [sp, #304]
 fmul	v1.4s, v1.4s, v0.4s
 ldr	q0, [sp, #368]
 fadd	v0.4s, v1.4s, v0.4s
 str	q0, [sp, #416]
 ldr	x0, [sp, #56]
 str	x0, [sp, #536]
 ldr	q0, [sp, #416]
 str	q0, [sp, #336]
 ldr	x0, [sp, #536]
 ldr	q0, [sp, #336]
 str	q0, [x0]
 nop
 add	sp, sp, #0x240
 ret
 nop
 .word	0x00000000
 .word	0x00000004
 .word	0x00000002
 .word	0x00000006
 .word	0x00000001
 .word	0x00000005
 .word	0x00000003
 .word	0x00000007
 .word	0x00000001
 .word	0x00000000
 .word	0x00000003
 .word	0x00000002
main:
 stp	x29, x30, [sp, #-32]!
 mov	x29, sp
 str	wzr, [x29, #28]
 ldr	w0, [x29, #28]
 cmp	w0, #0x3
 b.gt	400e00 <main+0xc0>
 ldr	w0, [x29, #28]
 add	w0, w0, #0x3
 scvtf	s0, w0
 adrp	x0, 412000 <__libc_start_main@GLIBC_2.17>
 add	x1, x0, #0x178
 ldrsw	x0, [x29, #28]
 lsl	x0, x0, #3
 add	x0, x1, x0
 str	s0, [x0]
 ldr	w0, [x29, #28]
 sub	w0, w0, #0x4
 scvtf	s0, w0
 adrp	x0, 412000 <__libc_start_main@GLIBC_2.17>
 add	x1, x0, #0x178
 ldrsw	x0, [x29, #28]
 lsl	x0, x0, #3
 add	x0, x1, x0
 add	x0, x0, #0x4
 str	s0, [x0]
 ldr	w0, [x29, #28]
 add	w0, w0, #0x4
 scvtf	s0, w0
 adrp	x0, 9ca8000 <A+0x9895e88>
 add	x1, x0, #0x9a0
 ldrsw	x0, [x29, #28]
 lsl	x0, x0, #3
 add	x0, x1, x0
 str	s0, [x0]
 ldr	w0, [x29, #28]
 add	w0, w0, #0x1
 scvtf	s0, w0
 adrp	x0, 9ca8000 <A+0x9895e88>
 add	x1, x0, #0x9a0
 ldrsw	x0, [x29, #28]
 lsl	x0, x0, #3
 add	x0, x1, x0
 add	x0, x0, #0x4
 str	s0, [x0]
 ldr	w0, [x29, #28]
 add	w0, w0, #0x1
 str	w0, [x29, #28]
 b	400d4c <main+0xc>
 adrp	x0, 412000 <__libc_start_main@GLIBC_2.17>
 add	x3, x0, #0x50
 adrp	x0, 1353f000 <B+0x9896660>
 add	x2, x0, #0x1c8
 adrp	x0, 9ca8000 <A+0x9895e88>
 add	x1, x0, #0x9a0
 adrp	x0, 412000 <__libc_start_main@GLIBC_2.17>
 add	x0, x0, #0x178
 bl	400880 <V_complex_mul(Complex*, Complex*, Complex*, float*)>
 str	wzr, [x29, #24]
 ldr	w0, [x29, #24]
 cmp	w0, #0x1
 b.gt	400ec0 <main+0x180>
 adrp	x0, 1353f000 <B+0x9896660>
 add	x1, x0, #0x1c8
 ldrsw	x0, [x29, #24]
 lsl	x0, x0, #3
 add	x0, x1, x0
 ldr	s0, [x0]
 adrp	x0, 412000 <__libc_start_main@GLIBC_2.17>
 add	x0, x0, #0x60
 bl	400730 <std::ostream::operator<<(float)@plt>
 mov	x2, x0
 adrp	x0, 400000 <_init-0x688>
 add	x0, x0, #0xff0
 mov	x1, x0
 mov	x0, x2
 bl	400710 <std::basic_ostream<char, std::char_traits<char> >& std::operator<< <std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*)@plt>
 mov	x2, x0
 adrp	x0, 1353f000 <B+0x9896660>
 add	x1, x0, #0x1c8
 ldrsw	x0, [x29, #24]
 lsl	x0, x0, #3
 add	x0, x1, x0
 add	x0, x0, #0x4
 ldr	s0, [x0]
 mov	x0, x2
 bl	400730 <std::ostream::operator<<(float)@plt>
 mov	x2, x0
 adrp	x0, 400000 <_init-0x688>
 add	x0, x0, #0xff0
 mov	x1, x0
 mov	x0, x2
 bl	400710 <std::basic_ostream<char, std::char_traits<char> >& std::operator<< <std::char_traits<char> >(std::basic_ostream<char, std::char_traits<char> >&, char const*)@plt>
 ldr	w0, [x29, #24]
 add	w0, w0, #0x1
 str	w0, [x29, #24]
 b	400e28 <main+0xe8>
 mov	w0, #0x0                   	// #0
 ldp	x29, x30, [sp], #32
 ret
int)::
 stp	x29, x30, [sp, #-32]!
 mov	x29, sp
 str	w0, [x29, #28]
 str	w1, [x29, #24]
 ldr	w0, [x29, #28]
 cmp	w0, #0x1
 b.ne	400f20 <__static_initialization_and_destruction_0(int, int)+0x54>  // b.any
 ldr	w1, [x29, #24]
 mov	w0, #0xffff                	// #65535
 cmp	w1, w0
 b.ne	400f20 <__static_initialization_and_destruction_0(int, int)+0x54>  // b.any
 adrp	x0, 1cdd5000 <C+0x9895e38>
 add	x0, x0, #0x9f0
 bl	4006e0 <std::ios_base::Init::Init()@plt>
 adrp	x0, 412000 <__libc_start_main@GLIBC_2.17>
 add	x2, x0, #0x48
 adrp	x0, 1cdd5000 <C+0x9895e38>
 add	x1, x0, #0x9f0
 adrp	x0, 400000 <_init-0x688>
 add	x0, x0, #0x700
 bl	4006f0 <__cxa_atexit@plt>
 nop
 ldp	x29, x30, [sp], #32
 ret
_GLOBAL__sub_I_A:
 stp	x29, x30, [sp, #-16]!
 mov	x29, sp
 mov	w1, #0xffff                	// #65535
 mov	w0, #0x1                   	// #1
 bl	400ecc <__static_initialization_and_destruction_0(int, int)>
 ldp	x29, x30, [sp], #16
 ret