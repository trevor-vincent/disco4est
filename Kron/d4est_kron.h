#ifndef D4EST_KRON_H
#define D4EST_KRON_H 

/* This file was automatically generated.  Do not edit! */
void d4est_kron_IoIoVEC_TRANSx_SQR(double *IoIoVEC_TRANSx,double *VEC,double *X,int N);
void d4est_kron_IoIoVECx_SQR(double *IoIoVECx,double *VEC,double *X,int N);
void d4est_kron_IoVEC_TRANSoIx_SQR(double *IoVEC_TRANSoIx,double *VEC,double *X,int N);
void d4est_kron_IoVECoIx_SQR(double *IoVECoIx,double *VEC,double *X,int N);
void d4est_kron_VEC_TRANSoIoIx_SQR(double *VEC_TRANSoIoIx,double *VEC,double *X,int N);
void d4est_kron_VECoIoIx_SQR(double *VECoIoIx,double *VEC,double *X,int N);
void d4est_kron_IoIoMATx_SQR(double *IoIoMATx,double *MAT,double *X,int N);
void d4est_kron_IoMAToIx_SQR(double *IoMAToIx,double *MAT,double *X,int N);
void d4est_kron_MAToIoIx_SQR(double *MAToIoIx,double *MAT,double *X,int N);
void d4est_kron_A1A2A3x_nonsqr(double *A1A2A3x,double *A1,double *A2,double *A3,double *X,int a1_rows,int a1_cols,int a2_rows,int a2_cols,int a3_rows,int a3_cols);
void d4est_kron_IoVEC_TRANSx_SQR(double *IoVEC_TRANSx,double *VEC,double *X,int N);
void d4est_kron_IoVECx_SQR(double *IoVECx,double *VEC,double *X,int N);
void d4est_kron_VEC_TRANSoIx_SQR(double *VEC_TRANSoIx,double *VEC,double *X,int N);
void d4est_kron_VECoIx_SQR(double *VECoIx,double *VEC,double *X,int N);
void d4est_kron_IoMATx_SQR(double *IoMATx,double *MAT,double *X,int N);
void d4est_kron_MAToIx_SQR(double *MAToIx,double *MAT,double *X,int N);
void d4est_kron_A1A2x_nonsqr(double *A1A2x,double *A1,double *A2,double *X,int a1_rows,int a1_cols,int a2_rows,int a2_cols);
void d4est_kron_AoBoC(double *A,double *B,double *C,double *D,int a_rows,int a_cols,int b_rows,int b_cols,int c_rows,int c_cols);
void d4est_kron_vec_o_vec_o_vec_dot_x(double *vec,double *x,int vec_size,double *vecvecvec_dot_x);
void d4est_kron_oneover_vec_o_vec_o_vec_dot_oneover_x_dot_y(double *vec,double *x,double *y,int vec_size,double *vecvecvec_dot_xy);
void d4est_kron_vec_o_vec_o_vec_dot_wxyz(double *vec,double *w,double *x,double *y,double *z,int vec_size,double *vecvecvec_dot_wxyz);
void d4est_kron_vec_o_vec_o_vec_dot_xy(double *vec,double *x,double *y,int vec_size,double *vecvecvec_dot_xy);
void d4est_kron_oneover_vec_dot_oneover_x_dot_y(double *vec,double *x,double *y,int vec_size,double *vec_dot_xy);
void d4est_kron_vec_dot_x(double *vec,double *x,int vec_size,double *vec_dot_x);
void d4est_kron_vec_dot_wxyz(double *vec,double *w,double *x,double *y,double *z,int vec_size,double *vec_dot_wxyz);
void d4est_kron_vec_dot_xy(double *vec,double *x,double *y,int vec_size,double *vec_dot_xy);
void d4est_kron_oneover_vec_o_vec_dot_oneover_x_dot_y(double *vec,double *x,double *y,int vec_size,double *vecvec_dot_xy);
void d4est_kron_vec_o_vec_dot_wxyz(double *vec,double *w,double *x,double *y,double *z,int vec_size,double *vecvec_dot_wxyz);
void d4est_kron_vec_o_vec_dot_xy(double *vec,double *x,double *y,int vec_size,double *vecvec_dot_xy);
void d4est_kron_vec1_o_vec2_o_vec3_dot_wxyz(double *vec1,double *vec2,double *vec3,double *w,double *x,double *y,double *z,int vec1_size,int vec2_size,int vec3_size,double *vec1vec2vec3_dot_x);
void d4est_kron_vec1_o_vec2_o_vec3_dot_xy(double *vec1,double *vec2,double *vec3,double *x,double *y,int vec1_size,int vec2_size,int vec3_size,double *vec1vec2vec3_dot_x);
void d4est_kron_vec1_o_vec2_o_vec3_dot_x(double *vec1,double *vec2,double *vec3,double *x,int vec1_size,int vec2_size,int vec3_size,double *vec1vec2vec3_dot_x);
void d4est_kron_vec1_o_vec2_dot_wxyz(double *vec1,double *vec2,double *w,double *x,double *y,double *z,int vec1_size,int vec2_size,double *vec1vec2_dot_x);
void d4est_kron_vec1_o_vec2_dot_xy(double *vec1,double *vec2,double *x,double *y,int vec1_size,int vec2_size,double *vec1vec2_dot_x);
void d4est_kron_vec1_o_vec2_dot_x(double *vec1,double *vec2,double *x,int vec1_size,int vec2_size,double *vec1vec2_dot_x);
void d4est_kron_vec_o_vec_dot_x(double *vec,double *x,int vec_size,double *vecvec_dot_x);
void d4est_kron_AoB(double *A,double *B,double *C,int a_rows,int a_cols,int b_rows,int b_cols);

#endif
