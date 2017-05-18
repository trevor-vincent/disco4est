#ifndef D4EST_GEOMETRY_CUBED_SPHERE_INVERSEMAP_QUADRATURE_H
#define D4EST_GEOMETRY_CUBED_SPHERE_INVERSEMAP_QUADRATURE_H 
#include <util.h>

/* jacobian goes as ~1/(c1+c2t)^4 */
typedef struct {
  long double cmin;
  long double cmax;
  long double R2;
  long double R1;
}inversemap_quadrature_params_t;

void
inversemap_aa_and_bb
(
 int n,
 long double* aa,
 long double* bb,
 void* user
)
{
  long double cmin = ((inversemap_quadrature_params_t*)user)->cmin;
  long double cmax = ((inversemap_quadrature_params_t*)user)->cmax;
  long double R1 = ((inversemap_quadrature_params_t*)user)->R1;
  long double R2 = ((inversemap_quadrature_params_t*)user)->R2;
  long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1);
  long double c2 = -((R2-R1)*(cmax-cmin));
  
  if (n == 1){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    bb[0] = 0;
  }
  else if (n == 2){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    aa[1] = (-9*powl(c1,5)*c2 + 5*c1*powl(c2,5) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*atanhl(c2/c1))/(3*powl(c1,2)*powl(c2,4) + powl(c2,6));
    bb[0] = 0;
    bb[1] = (3*powl(c1 - c2,2)*powl(c1 + c2,2))/powl(3*powl(c1,2) + powl(c2,2),2);
  }
  else if (n == 3){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    aa[1] = (-9*powl(c1,5)*c2 + 5*c1*powl(c2,5) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*atanhl(c2/c1))/(3*powl(c1,2)*powl(c2,4) + powl(c2,6));
    aa[2] = ((c1 - c2)*(c1 + c2)*(-9*c1*powl(c2,3)*(powl(c1,2) - 2*powl(c2,2)) + (27*powl(c1,4)*powl(c2,2) - 39*powl(c1,2)*powl(c2,4) - 2*powl(c2,6))*atanhl(c2/c1) + (-27*powl(c1,5)*c2 + 24*powl(c1,3)*powl(c2,3) + 7*c1*powl(c2,5))*powl(atanhl(c2/c1),2) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,4)*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2)));
    bb[0] = 0;
    bb[1] = (3*powl(c1 - c2,2)*powl(c1 + c2,2))/powl(3*powl(c1,2) + powl(c2,2),2);
    bb[2] = -(((c1 - c2)*(c1 + c2)*(3*powl(c1,2) + powl(c2,2))*(3*powl(c1,2)*powl(c2,2) - 4*powl(c2,4) + (-6*powl(c1,3)*c2 + 6*c1*powl(c2,3))*atanhl(c2/c1) + (3*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) - powl(c2,4))*powl(atanhl(c2/c1),2)))/powl(c2,8));
  }
  else if (n == 4){
    aa[0] = (-4*c1*c2)/(3*powl(c1,2) + powl(c2,2));
    aa[1] = (-9*powl(c1,5)*c2 + 5*c1*powl(c2,5) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*atanhl(c2/c1))/(3*powl(c1,2)*powl(c2,4) + powl(c2,6));
    aa[2] = ((c1 - c2)*(c1 + c2)*(-9*c1*powl(c2,3)*(powl(c1,2) - 2*powl(c2,2)) + (27*powl(c1,4)*powl(c2,2) - 39*powl(c1,2)*powl(c2,4) - 2*powl(c2,6))*atanhl(c2/c1) + (-27*powl(c1,5)*c2 + 24*powl(c1,3)*powl(c2,3) + 7*c1*powl(c2,5))*powl(atanhl(c2/c1),2) + (powl(c1,2) - powl(c2,2))*powl(3*powl(c1,2) + powl(c2,2),2)*powl(atanhl(c2/c1),3)))/(powl(c2,4)*(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2)));
    aa[3] = -(((c1 - c2)*(c1 + c2)*(-9*(3*powl(c1,5)*powl(c2,5) + 16*powl(c1,3)*powl(c2,7) - 16*c1*powl(c2,9)) + powl(c2,4)*(135*powl(c1,6) + 351*powl(c1,4)*powl(c2,2) - 576*powl(c1,2)*powl(c2,4) + 112*powl(c2,6))*atanhl(c2/c1) - 2*(135*powl(c1,7)*powl(c2,3) - 189*powl(c1,3)*powl(c2,7) + 56*c1*powl(c2,9))*powl(atanhl(c2/c1),2) + 2*powl(c2,2)*(135*powl(c1,8) - 315*powl(c1,6)*powl(c2,2) + 369*powl(c1,4)*powl(c2,4) - 245*powl(c1,2)*powl(c2,6) + 56*powl(c2,8))*powl(atanhl(c2/c1),3) + (-135*powl(c1,9)*c2 + 540*powl(c1,7)*powl(c2,3) - 882*powl(c1,5)*powl(c2,5) + 652*powl(c1,3)*powl(c2,7) - 175*c1*powl(c2,9))*powl(atanhl(c2/c1),4) + 27*powl(powl(c1,2) - powl(c2,2),4)*(powl(c1,2) + powl(c2,2))*powl(atanhl(c2/c1),5) - 36*c1*c2*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),6)))/((-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2))*(-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4))));
    bb[0] = 0;
    bb[1] = (3*powl(c1 - c2,2)*powl(c1 + c2,2))/powl(3*powl(c1,2) + powl(c2,2),2);
    bb[2] = -(((c1 - c2)*(c1 + c2)*(3*powl(c1,2) + powl(c2,2))*(3*powl(c1,2)*powl(c2,2) - 4*powl(c2,4) + (-6*powl(c1,3)*c2 + 6*c1*powl(c2,3))*atanhl(c2/c1) + (3*powl(c1,4) - 2*powl(c1,2)*powl(c2,2) - powl(c2,4))*powl(atanhl(c2/c1),2)))/powl(c2,8));
    bb[3] = (-27*powl(c1,4)*powl(c2,4) + 12*powl(c1,2)*powl(c2,6) + 16*powl(c2,8) + 24*c1*powl(c2,3)*(3*powl(c1,4) - 4*powl(c1,2)*powl(c2,2) + powl(c2,4))*atanhl(c2/c1) + (-54*powl(c1,6)*powl(c2,2) + 120*powl(c1,4)*powl(c2,4) - 94*powl(c1,2)*powl(c2,6) + 28*powl(c2,8))*powl(atanhl(c2/c1),2) + 9*powl(powl(c1,2) - powl(c2,2),4)*powl(atanhl(c2/c1),4))/(3.*powl(-3*powl(c1,2)*powl(c2,2) + 4*powl(c2,4) + 6*c1*c2*(powl(c1,2) - powl(c2,2))*atanhl(c2/c1) + (-3*powl(c1,4) + 2*powl(c1,2)*powl(c2,2) + powl(c2,4))*powl(atanhl(c2/c1),2),2));
  }
  /* else if (){ */


/* if (n = 5){ */

/* aa[0] = (-4*c1*c2)/(3*Power(c1,2) + Power(c2,2)); */

/* aa[1] = (-9*Power(c1,5)*c2 + 5*c1*Power(c2,5) + (Power(c1,2) - Power(c2,2))*Power(3*Power(c1,2) + Power(c2,2),2)*ArcTanh(c2/c1))/(3*Power(c1,2)*Power(c2,4) + Power(c2,6)); */

/* aa[2] = ((c1 - c2)*(c1 + c2)*(-9*c1*Power(c2,3)*(Power(c1,2) - 2*Power(c2,2)) + (27*Power(c1,4)*Power(c2,2) - 39*Power(c1,2)*Power(c2,4) - 2*Power(c2,6))*ArcTanh(c2/c1) + (-27*Power(c1,5)*c2 + 24*Power(c1,3)*Power(c2,3) + 7*c1*Power(c2,5))*Power(ArcTanh(c2/c1),2) + (Power(c1,2) - Power(c2,2))*Power(3*Power(c1,2) + Power(c2,2),2)*Power(ArcTanh(c2/c1),3)))/(Power(c2,4)*(-3*Power(c1,2)*Power(c2,2) + 4*Power(c2,4) + 6*c1*c2*(Power(c1,2) - Power(c2,2))*ArcTanh(c2/c1) + (-3*Power(c1,4) + 2*Power(c1,2)*Power(c2,2) + Power(c2,4))*Power(ArcTanh(c2/c1),2))); */

/* aa[3] = -(((c1 - c2)*(c1 + c2)*(-9*(3*Power(c1,5)*Power(c2,5) + 16*Power(c1,3)*Power(c2,7) - 16*c1*Power(c2,9)) + Power(c2,4)*(135*Power(c1,6) + 351*Power(c1,4)*Power(c2,2) - 576*Power(c1,2)*Power(c2,4) + 112*Power(c2,6))*ArcTanh(c2/c1) - 2*(135*Power(c1,7)*Power(c2,3) - 189*Power(c1,3)*Power(c2,7) + 56*c1*Power(c2,9))*Power(ArcTanh(c2/c1),2) + 2*Power(c2,2)*(135*Power(c1,8) - 315*Power(c1,6)*Power(c2,2) + 369*Power(c1,4)*Power(c2,4) - 245*Power(c1,2)*Power(c2,6) + 56*Power(c2,8))*Power(ArcTanh(c2/c1),3) + (-135*Power(c1,9)*c2 + 540*Power(c1,7)*Power(c2,3) - 882*Power(c1,5)*Power(c2,5) + 652*Power(c1,3)*Power(c2,7) - 175*c1*Power(c2,9))*Power(ArcTanh(c2/c1),4) + 27*Power(Power(c1,2) - Power(c2,2),4)*(Power(c1,2) + Power(c2,2))*Power(ArcTanh(c2/c1),5) - 36*c1*c2*Power(Power(c1,2) - Power(c2,2),4)*Power(ArcTanh(c2/c1),6)))/((-3*Power(c1,2)*Power(c2,2) + 4*Power(c2,4) + 6*c1*c2*(Power(c1,2) - Power(c2,2))*ArcTanh(c2/c1) + (-3*Power(c1,4) + 2*Power(c1,2)*Power(c2,2) + Power(c2,4))*Power(ArcTanh(c2/c1),2))*(-27*Power(c1,4)*Power(c2,4) + 12*Power(c1,2)*Power(c2,6) + 16*Power(c2,8) + 24*c1*Power(c2,3)*(3*Power(c1,4) - 4*Power(c1,2)*Power(c2,2) + Power(c2,4))*ArcTanh(c2/c1) + (-54*Power(c1,6)*Power(c2,2) + 120*Power(c1,4)*Power(c2,4) - 94*Power(c1,2)*Power(c2,6) + 28*Power(c2,8))*Power(ArcTanh(c2/c1),2) + 9*Power(Power(c1,2) - Power(c2,2),4)*Power(ArcTanh(c2/c1),4)))); */

/* aa[4] = -((81*Power(c2,2)*Power(-Power(c1,2) + Power(c2,2),5)*(45315*Power(c1,10) - 55485*Power(c1,8)*Power(c2,2) + 23902*Power(c1,6)*Power(c2,4) - 426*Power(c1,4)*Power(c2,6) + 111*Power(c1,2)*Power(c2,8) + 23*Power(c2,10))*Power(ArcTanh(c2/c1),7) + 306180*Power(c1,6)*Power(Power(c1,2) - Power(c2,2),8)*Power(ArcTanh(c2/c1),9) + 486*c1*Power(Power(c1,2) - Power(c2,2),8)*Power(ArcTanh(c2/c1),8)*(-10*Power(c1,2)*Power(c2,3) - 6*Power(c2,5) + 315*Power(c1,5)*Log(c1 - c2) - 315*Power(c1,5)*Log(c1 + c2)) + 9*Power(c2,3)*Power(Power(c1,2) - Power(c2,2),2)*Power(ArcTanh(c2/c1),5)*(1029105*Power(c1,14)*c2 - 2534220*Power(c1,12)*Power(c2,3) + 2920905*Power(c1,10)*Power(c2,5) - 1812408*Power(c1,8)*Power(c2,7) + 526595*Power(c1,6)*Power(c2,9) - 22092*Power(c1,4)*Power(c2,11) - 1869*Power(c1,2)*Power(c2,13) + 1504*Power(c2,15) + 90720*Power(c1,7)*Power(Power(c1,2) - Power(c2,2),3)*(3*Power(c1,2) - Power(c2,2))*Log(c1 - c2) - 90720*Power(c1,7)*Power(Power(c1,2) - Power(c2,2),3)*(3*Power(c1,2) - Power(c2,2))*Log(c1 + c2)) + 3*c1*Power(c2,8)*(-1215*Power(c1,12)*c2 - 2835*Power(c1,10)*Power(c2,3) + 35775*Power(c1,8)*Power(c2,5) - 94221*Power(c1,6)*Power(c2,7) + 111792*Power(c1,4)*Power(c2,9) - 62608*Power(c1,2)*Power(c2,11) + 13312*Power(c2,13) + 630*Power(c1,5)*Power(-27*Power(c1,4) + 12*Power(c1,2)*Power(c2,2) + 16*Power(c2,4),2)*Log(c1 - c2) - 630*Power(c1,5)*Power(-27*Power(c1,4) + 12*Power(c1,2)*Power(c2,2) + 16*Power(c2,4),2)*Log(c1 + c2)) - 27*c1*Power(c2,2)*Power(Power(c1,2) - Power(c2,2),5)*Power(ArcTanh(c2/c1),6)*(-180495*Power(c1,8)*c2 + 58500*Power(c1,6)*Power(c2,3) + 906*Power(c1,4)*Power(c2,5) + 1772*Power(c1,2)*Power(c2,7) - 1083*Power(c2,9) + 2520*(27*Power(c1,9) - 33*Power(c1,7)*Power(c2,2) + 14*Power(c1,5)*Power(c2,4))*Log(c1 - c2) - 2520*(27*Power(c1,9) - 33*Power(c1,7)*Power(c2,2) + 14*Power(c1,5)*Power(c2,4))*Log(c1 + c2)) + c1*Power(c2,6)*(Power(c1,2) - Power(c2,2))*Power(ArcTanh(c2/c1),2)*(-14773185*Power(c1,12)*c2 + 11484180*Power(c1,10)*Power(c2,3) + 7136262*Power(c1,8)*Power(c2,5) - 4352724*Power(c1,6)*Power(c2,7) + 1293147*Power(c1,4)*Power(c2,9) - 472704*Power(c1,2)*Power(c2,11) + 47888*Power(c2,13) + 7560*Power(c1,5)*(2025*Power(c1,8) - 3375*Power(c1,6)*Power(c2,2) + 1350*Power(c1,4)*Power(c2,4) + 216*Power(c1,2)*Power(c2,6) - 224*Power(c2,8))*Log(c1 - c2) - 7560*Power(c1,5)*(2025*Power(c1,8) - 3375*Power(c1,6)*Power(c2,2) + 1350*Power(c1,4)*Power(c2,4) + 216*Power(c1,2)*Power(c2,6) - 224*Power(c2,8))*Log(c1 + c2)) + Power(c2,7)*ArcTanh(c2/c1)*(2781135*Power(c1,14)*c2 - 2442150*Power(c1,12)*Power(c2,3) - 3187188*Power(c1,10)*Power(c2,5) + 2795094*Power(c1,8)*Power(c2,7) - 735831*Power(c1,6)*Power(c2,9) + 1074816*Power(c1,4)*Power(c2,11) - 307696*Power(c1,2)*Power(c2,13) + 25600*Power(c2,15) - 90720*Power(c1,7)*(81*Power(c1,8) - 144*Power(c1,6)*Power(c2,2) + 27*Power(c1,4)*Power(c2,4) + 52*Power(c1,2)*Power(c2,6) - 16*Power(c2,8))*Log(c1 - c2) + 90720*Power(c1,7)*(81*Power(c1,8) - 144*Power(c1,6)*Power(c2,2) + 27*Power(c1,4)*Power(c2,4) + 52*Power(c1,2)*Power(c2,6) - 16*Power(c2,8))*Log(c1 + c2)) - Power(c2,5)*(Power(c1,2) - Power(c2,2))*Power(ArcTanh(c2/c1),3)*(-30745575*Power(c1,14)*c2 + 51359265*Power(c1,12)*Power(c2,3) - 20519406*Power(c1,10)*Power(c2,5) - 3545478*Power(c1,8)*Power(c2,7) + 3445749*Power(c1,6)*Power(c2,9) + 299109*Power(c1,4)*Power(c2,11) - 205136*Power(c1,2)*Power(c2,13) + 32432*Power(c2,15) + 181440*Power(c1,7)*(81*Power(c1,8) - 207*Power(c1,6)*Power(c2,2) + 201*Power(c1,4)*Power(c2,4) - 89*Power(c1,2)*Power(c2,6) + 14*Power(c2,8))*Log(c1 - c2) - 181440*Power(c1,7)*(81*Power(c1,8) - 207*Power(c1,6)*Power(c2,2) + 201*Power(c1,4)*Power(c2,4) - 89*Power(c1,2)*Power(c2,6) + 14*Power(c2,8))*Log(c1 + c2)) + 3*c1*Power(c2,4)*Power(Power(c1,2) - Power(c2,2),2)*Power(ArcTanh(c2/c1),4)*(-9840285*Power(c1,12)*c2 + 15380685*Power(c1,10)*Power(c2,3) - 9305010*Power(c1,8)*Power(c2,5) + 2008170*Power(c1,6)*Power(c2,7) - 304065*Power(c1,4)*Power(c2,9) + 151897*Power(c1,2)*Power(c2,11) - 26880*Power(c2,13) + 1260*Power(c1,5)*(1215*Power(c1,8) - 2970*Power(c1,6)*Power(c2,2) + 3375*Power(c1,4)*Power(c2,4) - 2028*Power(c1,2)*Power(c2,6) + 536*Power(c2,8))*Log(c1 - c2) - 1260*Power(c1,5)*(1215*Power(c1,8) - 2970*Power(c1,6)*Power(c2,2) + 3375*Power(c1,4)*Power(c2,4) - 2028*Power(c1,2)*Power(c2,6) + 536*Power(c2,8))*Log(c1 + c2)))/(Power(c2,2)*(-27*Power(c1,4)*Power(c2,4) + 12*Power(c1,2)*Power(c2,6) + 16*Power(c2,8) + 24*c1*Power(c2,3)*(3*Power(c1,4) - 4*Power(c1,2)*Power(c2,2) + Power(c2,4))*ArcTanh(c2/c1) + (-54*Power(c1,6)*Power(c2,2) + 120*Power(c1,4)*Power(c2,4) - 94*Power(c1,2)*Power(c2,6) + 28*Power(c2,8))*Power(ArcTanh(c2/c1),2) + 9*Power(Power(c1,2) - Power(c2,2),4)*Power(ArcTanh(c2/c1),4))*(Power(c2,4)*(405*Power(c1,8) - 2970*Power(c1,6)*Power(c2,2) + 5157*Power(c1,4)*Power(c2,4) - 3612*Power(c1,2)*Power(c2,6) + 1024*Power(c2,8)) - 12*c1*Power(c2,3)*(135*Power(c1,8) - 810*Power(c1,6)*Power(c2,2) + 1431*Power(c1,4)*Power(c2,4) - 1018*Power(c1,2)*Power(c2,6) + 262*Power(c2,8))*ArcTanh(c2/c1) + 18*Power(c2,2)*(135*Power(c1,10) - 630*Power(c1,8)*Power(c2,2) + 1008*Power(c1,6)*Power(c2,4) - 640*Power(c1,4)*Power(c2,6) + 93*Power(c1,2)*Power(c2,8) + 34*Power(c2,10))*Power(ArcTanh(c2/c1),2) - 540*c1*c2*Power(Power(c1,2) - Power(c2,2),4)*(3*Power(c1,2) + 2*Power(c2,2))*Power(ArcTanh(c2/c1),3) + 81*Power(Power(c1,2) - Power(c2,2),4)*(5*Power(c1,4) + 10*Power(c1,2)*Power(c2,2) + Power(c2,4))*Power(ArcTanh(c2/c1),4)))); */

/* bb[0] = 0; */

/* bb[1] = (3*Power(c1 - c2,2)*Power(c1 + c2,2))/Power(3*Power(c1,2) + Power(c2,2),2); */

/* bb[2] = -(((c1 - c2)*(c1 + c2)*(3*Power(c1,2) + Power(c2,2))*(3*Power(c1,2)*Power(c2,2) - 4*Power(c2,4) + (-6*Power(c1,3)*c2 + 6*c1*Power(c2,3))*ArcTanh(c2/c1) + (3*Power(c1,4) - 2*Power(c1,2)*Power(c2,2) - Power(c2,4))*Power(ArcTanh(c2/c1),2)))/Power(c2,8)); */

/* bb[3] = (-27*Power(c1,4)*Power(c2,4) + 12*Power(c1,2)*Power(c2,6) + 16*Power(c2,8) + 24*c1*Power(c2,3)*(3*Power(c1,4) - 4*Power(c1,2)*Power(c2,2) + Power(c2,4))*ArcTanh(c2/c1) + (-54*Power(c1,6)*Power(c2,2) + 120*Power(c1,4)*Power(c2,4) - 94*Power(c1,2)*Power(c2,6) + 28*Power(c2,8))*Power(ArcTanh(c2/c1),2) + 9*Power(Power(c1,2) - Power(c2,2),4)*Power(ArcTanh(c2/c1),4))/(3.*Power(-3*Power(c1,2)*Power(c2,2) + 4*Power(c2,4) + 6*c1*c2*(Power(c1,2) - Power(c2,2))*ArcTanh(c2/c1) + (-3*Power(c1,4) + 2*Power(c1,2)*Power(c2,2) + Power(c2,4))*Power(ArcTanh(c2/c1),2),2)); */

/* bb[4] = ((-3*Power(c1,2)*Power(c2,2) + 4*Power(c2,4) + 6*c1*c2*(Power(c1,2) - Power(c2,2))*ArcTanh(c2/c1) + (-3*Power(c1,4) + 2*Power(c1,2)*Power(c2,2) + Power(c2,4))*Power(ArcTanh(c2/c1),2))*(Power(c2,4)*(405*Power(c1,8) - 2970*Power(c1,6)*Power(c2,2) + 5157*Power(c1,4)*Power(c2,4) - 3612*Power(c1,2)*Power(c2,6) + 1024*Power(c2,8)) - 12*c1*Power(c2,3)*(135*Power(c1,8) - 810*Power(c1,6)*Power(c2,2) + 1431*Power(c1,4)*Power(c2,4) - 1018*Power(c1,2)*Power(c2,6) + 262*Power(c2,8))*ArcTanh(c2/c1) + 18*Power(c2,2)*(135*Power(c1,10) - 630*Power(c1,8)*Power(c2,2) + 1008*Power(c1,6)*Power(c2,4) - 640*Power(c1,4)*Power(c2,6) + 93*Power(c1,2)*Power(c2,8) + 34*Power(c2,10))*Power(ArcTanh(c2/c1),2) - 540*c1*c2*Power(Power(c1,2) - Power(c2,2),4)*(3*Power(c1,2) + 2*Power(c2,2))*Power(ArcTanh(c2/c1),3) + 81*Power(Power(c1,2) - Power(c2,2),4)*(5*Power(c1,4) + 10*Power(c1,2)*Power(c2,2) + Power(c2,4))*Power(ArcTanh(c2/c1),4)))/(15.*Power(-27*Power(c1,4)*Power(c2,4) + 12*Power(c1,2)*Power(c2,6) + 16*Power(c2,8) + 24*c1*Power(c2,3)*(3*Power(c1,4) - 4*Power(c1,2)*Power(c2,2) + Power(c2,4))*ArcTanh(c2/c1) + (-54*Power(c1,6)*Power(c2,2) + 120*Power(c1,4)*Power(c2,4) - 94*Power(c1,2)*Power(c2,6) + 28*Power(c2,8))*Power(ArcTanh(c2/c1),2) + 9*Power(Power(c1,2) - Power(c2,2),4)*Power(ArcTanh(c2/c1),4),2)); */

/* } */
/*   } */
  else {
    mpi_abort("[D4EST_ERROR]: Do not support n >= 5 yet\n");
  }

}

long double
inversemap_weight_fcn(long double x, void* user){
  long double cmin = ((inversemap_quadrature_params_t*)user)->cmin;
  long double cmax = ((inversemap_quadrature_params_t*)user)->cmax;
  long double R1 = ((inversemap_quadrature_params_t*)user)->R1;
  long double R2 = ((inversemap_quadrature_params_t*)user)->R2;
  long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1);
  long double c2 = -((R2-R1)*(cmax-cmin));
  return 1.0l/powl(c1 + c2*x,4);
}

long double
inversemap_moment_fcn(int n, void* user)
{
  long double cmin = ((inversemap_quadrature_params_t*)user)->cmin;
  long double cmax = ((inversemap_quadrature_params_t*)user)->cmax;
  long double R1 = ((inversemap_quadrature_params_t*)user)->R1;
  long double R2 = ((inversemap_quadrature_params_t*)user)->R2;
  long double c1 = -((R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1);
  long double c2 = -((R2-R1)*(cmax-cmin));
  
  if (!(c1 + c2 > 0.0l && c1 - c2 > 0.0l && n <= 20)){
    printf(" c1 = (R2-R1)*(cmax+cmin) - 4.0l*R2 + 2.0l*R1 = %Le\n", c1);
    printf(" c2 = (R2-R1)*(cmax-cmin) = %Le\n", c2);
    printf(" R2 = %Le\n", R2);
    printf(" R1 = %Le\n", R1);
    printf(" cmax = %Le\n", cmax);
    printf(" cmin = %Le\n", cmin);
    printf(" n = %d\n", n);
    mpi_abort("[D4EST_ERROR]: condition c1 + c2 > 0 && c1 - c2 > 0 && n <= 20 failed\n");
  }
  
  if (n == 0) return (2*(3*powl(c1,2) + powl(c2,2)))/(3.*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 1) return (-8*c1*c2)/(3.*powl(c1 - c2,3)*powl(c1 + c2,3));
  else if (n == 2) return (2*(powl(c1,2) + 3*powl(c2,2)))/(3.*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 3) return (2*((-3*powl(c1,5)*c2 + 8*powl(c1,3)*powl(c2,3) - 9*c1*powl(c2,5))/powl(powl(c1,2) - powl(c2,2),3) + 3*atanhl(c2/c1)))/(3.*powl(c2,4));
  else if (n == 4) return (-24*powl(c1,6)*c2 + 64*powl(c1,4)*powl(c2,3) - 54*powl(c1,2)*powl(c2,5) + 6*powl(c2,7) + 24*c1*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1))/(3.*powl(c2,5)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 5) return (2*c1*(-30*powl(c1,6)*c2 + 80*powl(c1,4)*powl(c2,3) - 66*powl(c1,2)*powl(c2,5) + 12*powl(c2,7) + 30*c1*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(3.*powl(c2,6)*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 6) return (2*(-60*powl(c1,8)*c2 + 160*powl(c1,6)*powl(c2,3) - 132*powl(c1,4)*powl(c2,5) + 27*powl(c1,2)*powl(c2,7) + powl(c2,9) + 60*powl(powl(c1,3) - c1*powl(c2,2),3)*atanhl(c2/c1)))/(3.*powl(c2,7)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 7) return (c1*(-210*powl(c1,8)*c2 + 560*powl(c1,6)*powl(c2,3) - 462*powl(c1,4)*powl(c2,5) + 96*powl(c1,2)*powl(c2,7) + 8*powl(c2,9) + 210*powl(powl(c1,3) - c1*powl(c2,2),3)*atanhl(c2/c1)))/(3.*powl(c2,8)*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 8) return (840*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 2*(-840*powl(c1,10)*c2 + 2240*powl(c1,8)*powl(c2,3) - 1848*powl(c1,6)*powl(c2,5) + 384*powl(c1,4)*powl(c2,7) + 41*powl(c1,2)*powl(c2,9) + 3*powl(c2,11) + 420*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(15.*powl(c1 - c2,3)*powl(c2,9)*powl(c1 + c2,3));
  else if (n == 9) return (-4*c1*(-2*c2*(-315*powl(c1,10) + 840*powl(c1,8)*powl(c2,2) - 693*powl(c1,6)*powl(c2,4) + 144*powl(c1,4)*powl(c2,6) + 16*powl(c1,2)*powl(c2,8) + 3*powl(c2,10)) + 315*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 315*powl(c1,5)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(15.*powl(c1 - c2,3)*powl(c2,10)*powl(c1 + c2,3));
  else if (n == 10) return (2*(-2520*powl(c1,12)*c2 + 6720*powl(c1,10)*powl(c2,3) - 5544*powl(c1,8)*powl(c2,5) + 1152*powl(c1,6)*powl(c2,7) + 128*powl(c1,4)*powl(c2,9) + 33*powl(c1,2)*powl(c2,11) + 3*powl(c2,13) + 2520*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(21.*powl(c2,11)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 11) return (-2*c1*(-3465*powl(c1,12)*c2 + 9240*powl(c1,10)*powl(c2,3) - 7623*powl(c1,8)*powl(c2,5) + 1584*powl(c1,6)*powl(c2,7) + 176*powl(c1,4)*powl(c2,9) + 48*powl(c1,2)*powl(c2,11) + 12*powl(c2,13) + 3465*powl(c1,7)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(21.*powl(c2,12)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 12) return (2*(-((-13860*powl(c1,14)*c2 + 36960*powl(c1,12)*powl(c2,3) - 30492*powl(c1,10)*powl(c2,5) + 6336*powl(c1,8)*powl(c2,7) + 704*powl(c1,6)*powl(c2,9) + 192*powl(c1,4)*powl(c2,11) + 69*powl(c1,2)*powl(c2,13) + 7*powl(c2,15))/powl(powl(c1,2) - powl(c2,2),3)) - 13860*powl(c1,9)*atanhl(c2/c1)))/(63.*powl(c2,13));
  else if (n == 13) return (4*c1*((-45045*powl(c1,14)*c2 + 120120*powl(c1,12)*powl(c2,3) - 99099*powl(c1,10)*powl(c2,5) + 20592*powl(c1,8)*powl(c2,7) + 2288*powl(c1,6)*powl(c2,9) + 624*powl(c1,4)*powl(c2,11) + 240*powl(c1,2)*powl(c2,13) + 70*powl(c2,15))/powl(powl(c1,2) - powl(c2,2),3) + 45045*powl(c1,9)*atanhl(c2/c1)))/(315.*powl(c2,14));
  else if (n == 14) return (2*(-180180*powl(c1,16) + 480480*powl(c1,14)*powl(c2,2) - 396396*powl(c1,12)*powl(c2,4) + 82368*powl(c1,10)*powl(c2,6) + 9152*powl(c1,8)*powl(c2,8) + 2496*powl(c1,6)*powl(c2,10) + 960*powl(c1,4)*powl(c2,12) + 415*powl(c1,2)*powl(c2,14) + 45*powl(c2,16)))/(495.*powl(c2,14)*powl(-powl(c1,2) + powl(c2,2),3)) - (728*powl(c1,11)*atanhl(c2/c1))/powl(c2,15);
  else if (n == 15) return (-2*c1*(-45045*powl(c1,16)*c2 + 120120*powl(c1,14)*powl(c2,3) - 99099*powl(c1,12)*powl(c2,5) + 20592*powl(c1,10)*powl(c2,7) + 2288*powl(c1,8)*powl(c2,9) + 624*powl(c1,6)*powl(c2,11) + 240*powl(c1,4)*powl(c2,13) + 112*powl(c1,2)*powl(c2,15) + 36*powl(c2,17) + 45045*powl(c1,11)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(99.*powl(c2,16)*powl(-powl(c1,2) + powl(c2,2),3));
  else if (n == 16) return (2*(720720*powl(c1,18)*c2 - 1921920*powl(c1,16)*powl(c2,3) + 1585584*powl(c1,14)*powl(c2,5) - 329472*powl(c1,12)*powl(c2,7) - 36608*powl(c1,10)*powl(c2,9) - 9984*powl(c1,8)*powl(c2,11) - 3840*powl(c1,6)*powl(c2,13) - 1792*powl(c1,4)*powl(c2,15) - 873*powl(c1,2)*powl(c2,17) - 99*powl(c2,19) + 360360*powl(c1,13)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 360360*powl(c1,13)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(1287.*powl(c1 - c2,3)*powl(c2,17)*powl(c1 + c2,3));
  else if (n == 17) return (8*c1*(-1531530*powl(c1,18)*c2 + 4084080*powl(c1,16)*powl(c2,3) - 3369366*powl(c1,14)*powl(c2,5) + 700128*powl(c1,12)*powl(c2,7) + 77792*powl(c1,10)*powl(c2,9) + 21216*powl(c1,8)*powl(c2,11) + 8160*powl(c1,6)*powl(c2,13) + 3808*powl(c1,4)*powl(c2,15) + 2016*powl(c1,2)*powl(c2,17) + 693*powl(c2,19) + 1531530*powl(c1,13)*powl(powl(c1,2) - powl(c2,2),3)*atanhl(c2/c1)))/(9009.*powl(c2,18)*powl(powl(c1,2) - powl(c2,2),3));
  else if (n == 18) return (2*(12252240*powl(c1,20)*c2 - 32672640*powl(c1,18)*powl(c2,3) + 26954928*powl(c1,16)*powl(c2,5) - 5601024*powl(c1,14)*powl(c2,7) - 622336*powl(c1,12)*powl(c2,9) - 169728*powl(c1,10)*powl(c2,11) - 65280*powl(c1,8)*powl(c2,13) - 30464*powl(c1,6)*powl(c2,15) - 16128*powl(c1,4)*powl(c2,17) - 8547*powl(c1,2)*powl(c2,19) - 1001*powl(c2,21) + 6126120*powl(c1,15)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) - 6126120*powl(c1,15)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2)))/(15015.*powl(c1 - c2,3)*powl(c2,19)*powl(c1 + c2,3));
  else if (n == 19) return (2*c1*c2*(-14549535*powl(c1,20) + 38798760*powl(c1,18)*powl(c2,2) - 32008977*powl(c1,16)*powl(c2,4) + 6651216*powl(c1,14)*powl(c2,6) + 739024*powl(c1,12)*powl(c2,8) + 201552*powl(c1,10)*powl(c2,10) + 77520*powl(c1,8)*powl(c2,12) + 36176*powl(c1,6)*powl(c2,14) + 19152*powl(c1,4)*powl(c2,16) + 11088*powl(c1,2)*powl(c2,18) + 4004*powl(c2,20)) - 14549535*powl(c1,16)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 - c2) + 14549535*powl(c1,16)*powl(powl(c1,2) - powl(c2,2),3)*logl(c1 + c2))/(15015.*powl(c1 - c2,3)*powl(c2,20)*powl(c1 + c2,3));
  else if (n == 20) return (2*(-58198140*powl(c1,22) + 155195040*powl(c1,20)*powl(c2,2) - 128035908*powl(c1,18)*powl(c2,4) + 26604864*powl(c1,16)*powl(c2,6) + 2956096*powl(c1,14)*powl(c2,8) + 806208*powl(c1,12)*powl(c2,10) + 310080*powl(c1,10)*powl(c2,12) + 144704*powl(c1,8)*powl(c2,14) + 76608*powl(c1,6)*powl(c2,16) + 44352*powl(c1,4)*powl(c2,18) + 25025*powl(c1,2)*powl(c2,20) + 3003*powl(c2,22)))/(51051.*powl(c2,20)*powl(-powl(c1,2) + powl(c2,2),3)) - (2280*powl(c1,17)*atanhl(c2/c1))/powl(c2,21);
  else {
    mpi_abort("n must be <= 20");
    return NAN;
  }
}
 
#endif
