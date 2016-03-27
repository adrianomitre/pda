#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<time.h>
typedef struct l141 PDA_Parameters;struct l141{unsigned short
frame_len;unsigned int sample_rate;unsigned int bit_depth;unsigned
short max_submult;double dynamic_range;double min_expected_f0;double
max_expected_f0;double min_rel_prom;};const static unsigned short
DEFAULT_FRAME_LEN=1024;const static unsigned short DEFAULT_SAMPLE_RATE
=44100;const static unsigned short DEFAULT_BIT_DEPTH=16;const static
unsigned short DEFAULT_MAX_SUBMULTIPLE=20;const static double
DEFAULT_DYNAMIC_RANGE=50;const static double DEFAULT_MIN_REL_PROM=
0.97;extern void set_default_values(PDA_Parameters* );extern unsigned
char prepare_pda(const PDA_Parameters* );extern void
setup_dynamic_range(const double);extern double pda(const int[]);
extern double pda_prom(const int[],double* );extern double pda_wahm(
const int[],double* );extern double pda_prom_wahm(const int[],double*
,double* );extern void finalize_pda();
#define l65( l27) ( ( l27)/ floor( l27) <= ceil( l27)/( l27) ? floor( \
 l27) : ceil( l27) )
typedef struct l157{double lj;double ll;}ly;typedef enum{l30,l49}l111
;typedef struct l158{l111 l12;union l162{ly*l20;ly* *l23;}l13;}l70;
typedef struct l160{double l0;double l3;unsigned short l25;l70*lq;}
l63;typedef struct l161{double li;double lh;}l7;typedef enum{l29,l56}
l47;enum{l121,l87};const static unsigned char l148=255;const static
double l39=6.283185307179586476925286766559;const static double l84=
3.1415926535897932384626433832795;enum{l83=8,l91=24};enum{l77=128,l96
=16384};enum{l89=128,l95=192000};const static unsigned char l101=144;
const static unsigned short l45=16;const static double l54=0.5;const
static double l114=0.10034333188799373173791296490816;const static
double l88=1.02930223664349;const static double l76=0.971531941153606
;const static double l81=4186.0090448096;const static double l33=20;
extern void l109(const ly lk[],const unsigned short l36);extern void
l147(void);extern double l129(const int l1[]);extern double l102(void
);extern void l149(void);extern void l132(void);extern void l118(void
);extern void l105(void);extern ly*l145(void);static unsigned short*
l35;static double*l60, *l52, *l66, *l67;static double*l53;static
double*l62, *lx;static unsigned short l9;static l7*ls, *lt;static
double l128;static double l100,l153;static PDA_Parameters le;static
double l98;static double*l14, *l22, *l18;static double l99,l69;static
unsigned short lr;static ly l58;static ly*lo;static unsigned short lv
;static ly*lp;static double l71;static l63*lg;static unsigned short
l21,l5,lu;static l47 l48=l29;static unsigned char l82(const char*l140
,const char*l74){fprintf(stderr,"\x25\x73\x3a\x20\x25\x73\x2e\n",l140
,l74);return l87;}static unsigned char l16(const char*l74){return l82
("\x70\x72\x65\x70\x61\x72\x65\x5f\x70\x64\x61",l74);}static void*lw(
const size_t l42){void*lk;if((lk=malloc(l42))!=NULL)return lk;fprintf
(stderr,"\x45\x72\x72\x6f\x72\x20\x61\x6c\x6c\x6f\x63\x61\x74\x69\x6e"
"\x67\x20\x6d\x65\x6d\x6f\x72\x79\x2e\x20\x41\x62\x6f\x72\x74\x69\x6e"
"\x67\x2e\x2e\x2e\n");exit(l87);}static void l136(const unsigned short
lf){unsigned short ld,la,lc,lz=lf>>1,l75;l35=(unsigned short* )lw(lf*
sizeof(unsigned short));for(ld=0;ld<lf;ld++){for(la=1,lc=lz,l75=0;lc>
0;la=(la<<1),lc=(lc>>1))if(ld&lc)l75+=la;l35[ld]=l75;}}static void
l113(void){free(l35);}static void l103(const unsigned short lf){
unsigned short lm;l60=(double* )lw((lf+1) *sizeof(double));for(lm=2;
lm<=lf;lm<<=1)l60[lm]=sin(l39/lm);}static void l133(void){free(l60);}
static void l139(const unsigned short lf){unsigned short lm;l52=(
double* )lw((lf+1) *sizeof(double));for(lm=2;lm<=lf;lm<<=1)l52[lm]=
cos(l39/lm);}static void l116(void){free(l52);}static l47 l122(
unsigned long lf){if(!lf)return l29;while(lf>2){if(lf%2)return l29;lf
/=2;}return l56;}static void l108(const unsigned short lf){const
unsigned short lz=lf/2;unsigned short ld;l66=(double* )lw((lz) *
sizeof(double));for(ld=0;ld<lz;ld++)l66[ld]=sin(l39*ld/lf);}static
void l152(void){free(l66);}static void l154(const unsigned short lf){
const unsigned short lz=lf/2;unsigned short ld;l67=(double* )lw((lz) *
sizeof(double));for(ld=0;ld<lz;ld++)l67[ld]=cos(l39*ld/lf);}static
void l138(void){free(l67);}static void l112(const unsigned short lf){
l136(lf);l103(lf);l139(lf);}static void l117(void){l113();l133();l116
();}static void l134(){l154(le.frame_len);l108(le.frame_len);ls=(l7* )lw
(le.frame_len/2*sizeof(l7));lt=(l7* )lw(l9*sizeof(l7));l112(le.
frame_len/2);}static void l131(void){l138();l152();free(ls);free(lt);
l117();}static void l144(l7 ls[],const double lt[],const unsigned
short lf){const unsigned short lz=lf/2;unsigned short la,lc;for(la=lc
=0;la<lz;la++,lc++){ls[la].li=lt[lc];ls[la].lh=lt[++lc];}}static void
l150(l7 lt[],const l7 ls[],const unsigned short lf){const unsigned
short lz=lf/2;unsigned short ld;lt[0].li=ls[0].lh;lt[0].lh=0;for(ld=1
;ld<lz;ld++){const double l86=l67[ld],l46=l66[ld];const unsigned short
l50=lz-ld;lt[ld].li=0.5* (ls[ld].li* (1+l46)+ls[l50].li* (1-l46)+l86*
(ls[ld].lh+ls[l50].lh));lt[ld].lh=0.5* (ls[ld].lh* (1+l46)+ls[l50].lh
 * (l46-1)-l86* (ls[ld].li-ls[l50].li));}lt[lz].li=ls[0].li-ls[0].lh;
lt[lz].lh=0;}static void l127(const l7 l80[],l7 l4[],const unsigned
short int lf){register l7 l28,l2,l10,l34;register unsigned short int
ld,lm,lc;for(ld=0;ld<lf;ld++){l4[l35[ld]].li=l80[ld].li;l4[l35[ld]].
lh=l80[ld].lh;}for(lm=2;lm<=lf;lm=(lm<<1)){l28.li=l52[lm];l28.lh=l60[
lm];for(ld=0;ld<lf;ld+=lm){l2.li=1;l2.lh=0;for(lc=0;lc<(lm>>1);lc++){
l10.li=(l2.li*l4[ld+lc+(lm>>1)].li)-(l2.lh*l4[ld+lc+(lm>>1)].lh);l10.
lh=(l2.li*l4[ld+lc+(lm>>1)].lh)+(l2.lh*l4[ld+lc+(lm>>1)].li);l34.li=
l4[ld+lc].li;l34.lh=l4[ld+lc].lh;l4[ld+lc].li=l34.li+l10.li;l4[ld+lc]
.lh=l34.lh+l10.lh;l4[ld+lc+(lm>>1)].li=l34.li-l10.li;l4[ld+lc+(lm>>1)]
.lh=l34.lh-l10.lh;l10.li=(l2.li*l28.li)-(l2.lh*l28.lh);l10.lh=(l2.li*
l28.lh)+(l2.lh*l28.li);l2.li=l10.li;l2.lh=l10.lh;}}}}static void l130
(){const unsigned short l42=le.frame_len;unsigned short la;for(la=0;
la<l42;la++){l53[la]=0.5+0.5*cos(l84* (2*la-l42)/l42);}}static double
l32(const double l40){double ll;if(l40==0){return 0.25;}else{ll=sin(
l84*l40)/(l39*l40* (1-l40*l40));return ll*ll;}}static double l142(
const double l90,const double l79,const double l93){double l15;if(l90
>l93){l15=sqrt(l79/l90);return(2*l15-1)/(l15+1)-1;}else{l15=sqrt(l93/
l79);return(2*l15-1)/(l15+1);}}static void l115(){unsigned short la;
double l19,l68;lp=NULL;lv=0;for(la=1;la<l9-1;la++){if(lx[la]>=l100&&
lx[la]>=lx[la-1]&&lx[la]>lx[la+1]){l19=l142(lx[la-1],lx[la],lx[la+1]);
if(l19>=-l54&&l19<l54){l68=10*log10(lx[la]/l32(l19))-l98;assert(l32(
l19)<=l32(0));assert(l32(l19)>=l32(l54));if(l68>0){++lv;lo[(lv)].ll=
l68;lo[(lv)].lj=la+l19;if(!lp||lo[(lv)].ll>lp->ll){lp=&lo[(lv)];}}}la
++;}}if(lv)l71=lo[(lv)].lj;}static void l155(const int l123[]){
unsigned short la;for(la=0;la<le.frame_len;la++){l62[la]=l123[la] *
l53[la];}}static void l126(void){unsigned short la;l144(lt,l62,le.
frame_len);l127(lt,ls,le.frame_len/2);l150(lt,ls,le.frame_len);for(la
=0;la<l9;la++){lx[la]=lt[la].li*lt[la].li+lt[la].lh*lt[la].lh;}lx[l9-
1]/=4;}static void l137(){unsigned short lb;double l15,l97;l14=(
double* )lw((lr+1) *sizeof(double));for(lb=1;lb<=((4)<=(lr)?(4):(lr));
lb++){l14[(lb)]=1;}for(lb=5;lb<=lr;lb++){l15=sqrt(lb/(lb-1)) * (lb-1);
l97=sqrt((lb+1)/lb) *lb;l14[(lb)]=log10(l97/l15)/l114;}}static void
l146(){register unsigned short lb;l22=(double* )lw((lr+1) *sizeof(
double));l18=(double* )lw((lr+1) *sizeof(double));for(lb=1;lb<=((16)<=
(lr)?(16):(lr));++lb){l22[(lb)]=lb*l76;l18[(lb)]=lb*l88;}if(lr>=17){
l22[(17)]=17*l76;l18[(17)]=17*sqrt((lb+1.0)/lb);}for(lb=18;lb<=lr;++
lb){l22[(lb)]=lb*sqrt((lb-1.0)/lb);l18[(lb)]=lb*sqrt((lb+1.0)/lb);}}
static void l120(){unsigned short la,lc,ld;l21=((floor(le.sample_rate
/(2*le.min_expected_f0)))<=(le.max_submult)?(floor(le.sample_rate/(2*
le.min_expected_f0))):(le.max_submult));lg=(l63* )lw((l21+1) *sizeof(
l63));for(la=1;la<=l21;++la){lg[(la)].lq=(l70* )lw((lr+1) *sizeof(l70
));for(lc=1;lc<=lr;++lc){lg[(la)].lq[(lc)].l12=l30;}}for(la=1;la<=((
l21)<=(le.max_submult/2)?(l21):(le.max_submult/2));++la){for(lc=1;lc
<=l45/2;++lc){for(ld=2;ld*la<=l21&&ld*lc<=((lr)<=(l45)?(lr):(l45));++
ld){if(lg[(la*ld)].lq[(lc*ld)].l12!=l49){lg[(la*ld)].lq[(lc*ld)].l12=
l49;lg[(la*ld)].lq[(lc*ld)].l13.l23=&lg[(la)].lq[(lc)].l13.l20;}}}}}
static void l110(){unsigned short la;for(la=1;la<=l21;++la){free(lg[(
la)].lq);}free(lg);}static ly*l31(const unsigned short ln,const
unsigned short lb){if(lg[(ln)].lq[(lb)].l12==l49){return* (lg[(ln)].
lq[(lb)].l13.l23);}else{return lg[(ln)].lq[(lb)].l13.l20;}}static void
l135(){unsigned short la,lc;double l85;if(lv==0){lu=0;return;}assert(
lp!=NULL);l5=ceil(lp->lj/l99);lu=((le.max_submult)<=(floor(lp->lj/l69
))?(le.max_submult):(floor(lp->lj/l69)));if(lu<l5){lu=0;return;}
assert(l5>0);assert(lu<=lr);assert(l5<=lu);for(la=1;la<=lu;++la){l85=
lp->lj/la;lg[(la)].l25=l65(l71/l85);for(lc=1;lc<=lg[(la)].l25;++lc){
if(lg[(la)].lq[(lc)].l12==l30){lg[(la)].lq[(lc)].l13.l20=&l58;}assert
(lg[(la)].lq[(lc)].l12==l30||l31(la,lc)->ll==0);assert(lg[(la)].lq[(
lc)].l12==l30||l31(la,lc)->lj==0.123456789);}}}static double l94(
const double l57,const double l23){return((l57)>=(l23)?(l57):(l23))/(
(l57)<=(l23)?(l57):(l23));}static void l72(const unsigned short ln,
const unsigned short lk,const unsigned short lb){if(lo[(lk)].ll>lg[(
ln)].lq[(lb)].l13.l20->ll||(lo[(lk)].ll==lg[(ln)].lq[(lb)].l13.l20->
ll&&l94(lo[(lk)].lj,lp->lj/ln*lb)<l94(lg[(ln)].lq[(lb)].l13.l20->lj,
lp->lj/ln*lb))){lg[(ln)].lq[(lb)].l13.l20=&lo[(lk)];}}static l47 l64(
const double l6,const unsigned short lb){if(l6>l22[(lb)]&&l6<l18[(lb)]
){return l56;}else{return l29;}}static void l107(){unsigned short la,
lc,lb;double l6,l37;for(la=1;la<=lu;++la){l37=(double)la/lp->lj;for(
lc=1;lc<=lv;++lc){assert(lo[(lc)].ll>0);l6=lo[(lc)].lj*l37;lb=l65(l6);
assert(lb<=lr);if(lg[(la)].lq[(lb)].l12==l30){if(l64(l6,lb)){l72(la,
lc,lb);}}}}}static void l159(){unsigned short la,lc,ld,lb;double l6,
l37;for(la=l5;la<=lu;++la){l37=(double)la/lp->lj;for(lc=lv;lc>0;--lc){
assert(lo[(lc)].ll>0);l6=lo[(lc)].lj*l37;lb=l65(l6);assert(lb<=lr);if
(lg[(la)].lq[(lb)].l12==l30){if(l64(l6,lb)){l72(la,lc,lb);for(ld=(l45
/lb)+1;ld*la<=lu;++ld){if(l64(ld*l6,ld*lb)){l72(ld*la,lc,ld*lb);}}}}}
}}static void l106(){unsigned short la,lc;double l43,l8,l59;for(la=l5
;la<=lu;++la){l43=l8=l59=0;for(lc=1;lc<=lg[(la)].l25;++lc){if(l31(la,
lc)->ll>0){l43+=l31(la,lc)->ll*l14[(lc)];l8+=l14[(lc)]+l59;l59=0;}
else{l59+=l14[(lc)];}}lg[(la)].l0=l43;lg[(la)].l3=l8>0?l43/l8:0;}}
static double l151(const unsigned short ln){double l8=0,l61=0,l51;
unsigned short lb;ly*l17;assert(ln>=l5);assert(ln<=lu);assert(lv>0);
assert(lp!=NULL);assert(ln<=lu);assert(lu>0);for(lb=1;lb<=lg[(ln)].
l25;++lb){l17=l31(ln,lb);if(l17==&l58){continue;}assert(l17->ll!=0.0
&&l17->lj!=0.123456789);assert(l17!=NULL);l51=l17->ll*l14[(lb)];
assert(l51>=0);l61+=(l17->lj/lb) *l51;assert(l64(l17->lj/(lp->lj/ln),
lb));l8+=l51;}assert(l8!=0);assert(l61/l8>(lp->lj/ln) *l76);assert(
l61/l8<(lp->lj/ln) *l88);return l61/l8;}void l147(){unsigned short la
;for(la=0;la<=l9;la++){fprintf(stderr,"\x70\x6f\x77\x5f\x73\x70\x65"
"\x63\x5b\x25\x64\x5d\x20\x3d\x20\x25\x67\n",la,lx[la]);}}double l129
(const int l1[]){double l38=0;unsigned short la;for(la=0;la<le.
frame_len;la++){l38+=(l1[la] *l1[la]);}return l38;}double l102(){
double l38=0;unsigned short la;for(la=0;la<=l9;la++){l38+=lx[la];}
return l38;}static void l143(unsigned lf){fprintf(stderr,"\x69\x6e"
"\x68\x61\x72\x6d\x28\x25\x64\x29\x3d\x20\x5b\x25\x67\x2c\x20\x25\x67"
"\x5d\n",lf,l22[(lf)],l18[(lf)]);}void l149(){unsigned short la;for(
la=1;la<=lr;++la){l143(la);}}static void l125(unsigned lf){fprintf(
stderr,"\x77\x28\x25\x64\x29\x3d\x20\x25\x67\n",lf,l14[lf-1]);}void
l132(){unsigned short la;for(la=1;la<=lr;++la){l125(la);}}ly*l145(){
return lp;}static void l104(){unsigned short la;lp=NULL;for(la=1;la<=
lv;++la){if(!lp||lo[(la)].ll>lp->ll){lp=&lo[(la)];}}}static void l119
(unsigned short la){if(lp==&lo[la]){fprintf(stderr,"\x70\x28\x25\x64"
"\x29\x3d\x20\x5b\x66\x3a\x20\x25\x67\x2c\x20\x25\x67\x5d\x20\x3c\x3d"
"\x3d\x3d\x20\x73\x74\x72\x6f\x6e\x67\x65\x73\x74\n",la+1,lo[la].lj,
lo[la].ll);}else{fprintf(stderr,"\x70\x28\x25\x64\x29\x3d\x20\x5b\x66"
"\x3a\x20\x25\x67\x2c\x20\x25\x67\x5d\n",la+1,lo[la].lj,lo[la].ll);}}
void l118(void){unsigned short la;for(la=0;la<lv;la++){l119(la);}}
static l47 l156(const unsigned short ln,const unsigned short lb){
return lg[(ln)].lq[(lb)].l12==l49?l56:l29;}static void l124(unsigned
short lf){const l63*ln=&lg[(lf)];ly*lk;unsigned short lb;fprintf(
stderr,"\x63\x61\x6e\x64\x28\x25\x64\x29\x20\x3d\x3e\x20\x67\x72\x73"
"\x5f\x66\x72\x71\x3a\x20\x25\x67\x2c\x20\x70\x72\x6f\x6d\x3a\x20\x25"
"\x67\x2c\x20\x77\x61\x68\x6d\x3a\x20\x25\x67\x20\x28\x68\x69\x67\x68"
"\x2e\x70\x6f\x73\x2e\x68\x2e\x25\x64\x29\n",lf,(lp->lj)/lf,ln->l0,ln
->l3,ln->l25);for(lb=1;lb<=lg[(lf)].l25;++lb){lk=l31(lf,lb);assert(lk
!=NULL);if(lk->ll!=0){fprintf(stderr,"\t\x68\x28\x25\x64\x29\x20\x3d"
"\x20\x5b\x25\x67\x2c\x20\x25\x67\x5d\x20\x25\x73\n",lb,lk->lj,lk->ll
,(l156(lf,lb)?"\x28\x72\x65\x66\x29":""));}}}void l105(){unsigned
short la;for(la=l5;la<=lu;++la){l124(la);}}void l109(const ly lk[],
const unsigned short l36){unsigned short la;assert(l36>0);assert(l36
<=l9/2);assert(lk!=NULL);for(la=1;la<=l36;++la){lo[(la)].lj=lk[la-1].
lj;lo[(la)].ll=lk[la-1].ll;}lv=l36;if(lv){l71=lo[(lv)].lj;l104();}}
void set_default_values(PDA_Parameters*lk){lk->dynamic_range=
DEFAULT_DYNAMIC_RANGE;lk->min_rel_prom=DEFAULT_MIN_REL_PROM;lk->
max_submult=DEFAULT_MAX_SUBMULTIPLE;lk->frame_len=DEFAULT_FRAME_LEN;
lk->sample_rate=DEFAULT_SAMPLE_RATE;lk->bit_depth=DEFAULT_BIT_DEPTH;
lk->min_expected_f0=0;lk->max_expected_f0=0;}void setup_dynamic_range
(const double l44){if(l44<=0||l44>l101){l82("\x70\x64\x61\x5f\x64\x79"
"\x6e\x61\x6d\x69\x63\x5f\x72\x61\x6e\x67\x65","\x64\x79\x6e\x61\x6d"
"\x69\x63\x20\x72\x61\x6e\x67\x65\x20\x6d\x75\x73\x74\x20\x73\x61\x74"
"\x69\x73\x66\x79\x20\x30\x20\x3c\x20\x44\x52\x20\x3c\x3d\x20\x31\x34"
"\x34\x20\x64\x42");}double l92=le.frame_len*le.frame_len*pow(4,le.
bit_depth-2);l98=10*log10(l92)-l44;l153=l92*pow(10,-0.1*l44);}
unsigned char prepare_pda(const PDA_Parameters*l73){const double l55=
l73->sample_rate/2.0;char l11[l148];double l26;if(l48){return l16(""
"\x63\x61\x6e\x6e\x6f\x74\x20\x69\x6e\x69\x74\x69\x61\x6c\x69\x7a\x65"
"\x20\x6d\x75\x6c\x74\x69\x70\x6c\x65\x20\x74\x69\x6d\x65\x73\x20\x63"
"\x6f\x6e\x73\x65\x63\x75\x74\x69\x76\x65\x6c\x79");}assert(l73!=NULL
);memcpy(&le,l73,sizeof(PDA_Parameters));if(le.sample_rate<l89||le.
sample_rate>l95){sprintf(l11,"\x73\x61\x6d\x70\x6c\x65\x5f\x72\x61"
"\x74\x65\x20\x6d\x75\x73\x74\x20\x62\x65\x20\x62\x65\x74\x77\x65\x65"
"\x6e\x20\x25\x64\x20\x61\x6e\x64\x20\x25\x64",l89,l95);return l16(
l11);}if(le.frame_len<l77||le.frame_len>l96||!l122(le.frame_len)){
sprintf(l11,"\x66\x72\x61\x6d\x65\x5f\x6c\x65\x6e\x20\x6d\x75\x73\x74"
"\x20\x62\x65\x20\x61\x20\x70\x6f\x77\x65\x72\x20\x6f\x66\x20\x32\x20"
"\x77\x69\x74\x68\x69\x6e\x20\x25\x64\x20\x61\x6e\x64\x20\x25\x64",
l77,l96);return l16(l11);}l26=2*le.sample_rate/(double)le.frame_len;
if(!le.min_expected_f0){le.min_expected_f0=l26<l33?l33:l26;}if(!le.
max_expected_f0){le.max_expected_f0=l55>l81?l81:l55;}if(le.
min_expected_f0<l33||le.min_expected_f0>l55){sprintf(l11,"\x6d\x69"
"\x6e\x5f\x65\x78\x70\x65\x63\x74\x65\x64\x5f\x66\x30\x20\x6d\x75\x73"
"\x74\x20\x6c\x69\x65\x20\x62\x65\x74\x77\x65\x65\x6e\x20\x25\x67\x20"
"\x48\x7a\x20\x61\x6e\x64\x20\x4e\x79\x71\x75\x69\x73\x74\x20\x66\x72"
"\x65\x71\x75\x65\x6e\x63\x79",l33);return l16(l11);}if(le.
max_expected_f0<l33||le.max_expected_f0>l55){sprintf(l11,"\x6d\x61"
"\x78\x5f\x65\x78\x70\x65\x63\x74\x65\x64\x5f\x66\x30\x20\x6d\x75\x73"
"\x74\x20\x6c\x69\x65\x20\x62\x65\x74\x77\x65\x65\x6e\x20\x25\x67\x20"
"\x48\x7a\x20\x61\x6e\x64\x20\x4e\x79\x71\x75\x69\x73\x74\x20\x66\x72"
"\x65\x71\x75\x65\x6e\x63\x79",l33);return l16(l11);}if(le.
min_expected_f0>=le.max_expected_f0){return l16("\x6d\x61\x78\x5f\x65"
"\x78\x70\x65\x63\x74\x65\x64\x5f\x66\x30\x20\x6d\x75\x73\x74\x20\x62"
"\x65\x20\x67\x72\x65\x61\x74\x65\x72\x20\x74\x68\x61\x6e\x20\x6d\x69"
"\x6e\x5f\x65\x78\x70\x65\x63\x74\x65\x64\x5f\x66\x30");}if(le.
min_expected_f0<l26){return l16("\x6d\x69\x6e\x5f\x65\x78\x70\x65\x63"
"\x74\x65\x64\x5f\x66\x30\x20\x63\x61\x6e\x6e\x6f\x74\x20\x62\x65\x20"
"\x73\x6d\x61\x6c\x6c\x65\x72\x20\x74\x68\x61\x6e\x20\x74\x68\x65\x20"
"\x66\x72\x65\x71\x75\x65\x6e\x63\x79\x20\x72\x65\x73\x6f\x6c\x75\x74"
"\x69\x6f\x6e\x2c\n" "\x77\x68\x69\x63\x68\x20\x69\x73\x20\x73\x61"
"\x6d\x70\x6c\x65\x20\x72\x61\x74\x65\x20\x64\x69\x76\x69\x64\x65\x64"
"\x20\x62\x79\x20\x68\x61\x6c\x66\x20\x74\x68\x65\x20\x66\x72\x61\x6d"
"\x65\x20\x6c\x65\x6e\x67\x74\x68");}if(le.dynamic_range<=0){return
l16("\x69\x6e\x76\x61\x6c\x69\x64\x20\x64\x65\x66\x5f\x6e\x6f\x72\x6d"
"\x5f\x66\x61\x63\x74\x6f\x72");}if(le.bit_depth%8||le.bit_depth<l83
||le.bit_depth>l91){sprintf(l11,"\x62\x69\x74\x5f\x64\x65\x70\x74\x68"
"\x20\x6d\x75\x73\x74\x20\x62\x65\x20\x61\x20\x6d\x75\x6c\x74\x69\x70"
"\x6c\x65\x20\x6f\x66\x20\x77\x65\x69\x67\x68\x74\x20\x77\x69\x74\x68"
"\x69\x6e\x20\x25\x64\x20\x61\x6e\x64\x20\x25\x64",l83,l91);return l16
(l11);}l69=2*le.min_expected_f0/l26;l99=2*le.max_expected_f0/l26;l53=
(double* )lw(le.frame_len*sizeof(double));l130();l62=(double* )lw(le.
frame_len*sizeof(double));l9=le.frame_len/2+1;lx=(double* )lw(l9*
sizeof(double));lo=(ly* )lw((l9/2+1) *sizeof(ly));setup_dynamic_range
(le.dynamic_range);l134();l128=l32(l54);lr=l65(le.sample_rate/(2*le.
min_expected_f0));l137();l146();l120();l58.lj=0.123456789;l58.ll=0.0;
l48=l56;return l121;}void finalize_pda(){if(l48){l131();free(lo);free
(lx);free(l62);free(l53);free(l18);free(l22);free(l14);l110();l48=l29
;}}double pda(const int l1[]){return pda_prom(l1,NULL);}double
pda_prom(const int l1[],double*l0){return pda_prom_wahm(l1,l0,NULL);}
double pda_wahm(const int l1[],double*l3){return pda_prom_wahm(l1,
NULL,l3);}double pda_prom_wahm(const int l1[],double*l0,double*l3){
unsigned short la,l41,l24;double l78;assert(l1!=NULL);l155(l1);l126();
l115();l135();if(lu==0){if(l0) *l0=0;if(l3) *l3=0;return 0;}l107();
l106();l41=l5;for(la=l5+1;la<=lu;++la){if(lg[(la)].l0>lg[(l41)].l0){
l41=la;}}l78=le.min_rel_prom*lg[(l41)].l0;l24=l41;for(la=l5;la<=lu;++
la){if(lg[(la)].l0>=l78&&lg[(la)].l3>lg[(l24)].l3){l24=la;}}if(l0) *
l0=lg[(l24)].l0;if(l3) *l3=lg[(l24)].l3;return l151(l24) *le.
sample_rate/le.frame_len;}
