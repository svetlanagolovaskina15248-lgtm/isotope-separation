// ЗАПУСК (Мак, main.c находится в загрузках):
// cd ~/Downloads

// clang -std=c11 -O2 main.c -o isotopes -lm

// ./isotopes

// ИЗМЕНЕНИЕ
// ./isotopes --B 0.02 --v 1500 --angle 15

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265358979323846
#define AMU 1.66053906660e-27
#define Q 1.6e-19
#define MAXPTS 6001

typedef struct { double x,y,z,vx,vy,vz; } State;
typedef struct { double mass,B,Ex,Ey,Ez,v,vz,angle; } Params;
typedef struct { double t,x,y,z; } Point;
typedef struct { Point p[MAXPTS]; int n; double hit,energy,err; int crossed; } Result;
typedef struct { double x,y; } XY;

enum { RK4, BORIS, EULER };
static const char *palette[]={"#1769aa","#168b72","#c48b19","#d34b48","#8053ad","#7a6e5d"};

static State add(State a,State b,double h) {
    State c={a.x+h*b.x,a.y+h*b.y,a.z+h*b.z,a.vx+h*b.vx,a.vy+h*b.vy,a.vz+h*b.vz}; return c;
}
static State rhs(State s,Params p) {
    double k=Q/(p.mass*AMU);
    State d={s.vx,s.vy,s.vz,k*(p.Ex+s.vy*p.B),k*(p.Ey-s.vx*p.B),k*p.Ez}; return d;
}
static State step_rk4(State s,Params p,double h) {
    State a=rhs(s,p),b=rhs(add(s,a,h/2),p),c=rhs(add(s,b,h/2),p),d=rhs(add(s,c,h),p);
    State t=add(s,a,h/6);t=add(t,b,h/3);t=add(t,c,h/3);return add(t,d,h/6);
}
static State step_euler(State s,Params p,double h) {return add(s,rhs(s,p),h);}
static State step_boris(State s,Params p,double h) {
    double k=Q/(p.mass*AMU),t=k*p.B*h/2,ss=2*t/(1+t*t);
    double xm=s.vx+k*p.Ex*h/2,ym=s.vy+k*p.Ey*h/2,zm=s.vz+k*p.Ez*h/2;
    double xp=xm+ym*t,yp=ym-xm*t;
    double xn=xm+yp*ss,yn=ym-xp*ss;
    double vx=xn+k*p.Ex*h/2,vy=yn+k*p.Ey*h/2,vz=zm+k*p.Ez*h/2;
    State o={s.x+(s.vx+vx)*h/2,s.y+(s.vy+vy)*h/2,s.z+(s.vz+vz)*h/2,vx,vy,vz};return o;
}
static double speed2(State s){return s.vx*s.vx+s.vy*s.vy+s.vz*s.vz;}
static State initial(Params p){
    double a=p.angle*PI/180;
    State s={0,0,0,p.v*sin(a),p.v*cos(a),p.vz};return s;
}
static State exact(Params p,double t){
    State s=initial(p);double w=Q*p.B/(p.mass*AMU),a=w*t;
    double vx=s.vx,vy=s.vy;
    s.vx=vx*cos(a)+vy*sin(a);s.vy=vy*cos(a)-vx*sin(a);
    s.x=(vx*sin(a)+vy*(1-cos(a)))/w;
    s.y=(vy*sin(a)-vx*(1-cos(a)))/w;
    s.z=s.vz*t;return s;
}
static int zero_e(Params p){return p.Ex==0 && p.Ey==0 && p.Ez==0;}
static Result *simulate(Params p,int scheme,int steps,double periods,const char *csv){
    if(steps<8||steps*periods>=MAXPTS) return NULL;
    Result *r=calloc(1,sizeof(*r));if(!r)return NULL;
    FILE *f=NULL;if(csv){f=fopen(csv,"w");if(!f){free(r);return NULL;}fprintf(f,"t_s,x_mm,y_mm,z_mm,vx_m_s,vy_m_s,vz_m_s\n");}
    double T=2*PI*p.mass*AMU/(Q*p.B),h=T/steps;
    int count=(int)round(steps*periods);
    State s=initial(p);double e0=speed2(s);r->hit=NAN;
    for(int i=0;i<=count;i++){
        double t=i*h;
        if(i<MAXPTS){r->p[i]=(Point){t,s.x*1e3,s.y*1e3,s.z*1e3};r->n++;}
        if(f)fprintf(f,"%.12g,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n",t,s.x*1e3,s.y*1e3,s.z*1e3,s.vx,s.vy,s.vz);
        if(zero_e(p)){
            State a=exact(p,t);
            double err=hypot(hypot(s.x-a.x,s.y-a.y),s.z-a.z)*1e3;
            if(err>r->err)r->err=err;
            double en=fabs(speed2(s)/e0-1);
            if(en>r->energy)r->energy=en;
        }
        if(i==count)break;
        State next=scheme==RK4?step_rk4(s,p,h):scheme==BORIS?step_boris(s,p,h):step_euler(s,p,h);
        if(!r->crossed && i>0 && s.y>0 && next.y<=0){
            double u=s.y/(s.y-next.y);
            r->hit=(s.x+u*(next.x-s.x))*1e3;r->crossed=1;
        }
        s=next;
    }
    if(f)fclose(f);return r;
}
static void svg_begin(FILE *f,const char *title,const char *subtitle){
    fprintf(f,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"1100\" height=\"740\" viewBox=\"0 0 1100 740\">\n");
    fprintf(f,"<rect width=\"1100\" height=\"740\" fill=\"white\"/>\n<g font-family=\"Arial, sans-serif\" fill=\"#172b4d\">\n");
    fprintf(f,"<text x=\"75\" y=\"48\" font-size=\"26\" font-weight=\"bold\">%s</text>\n",title);
    fprintf(f,"<text x=\"75\" y=\"80\" font-size=\"16\">%s</text>\n",subtitle);
}
static void svg_end(FILE *f){fprintf(f,"</g>\n</svg>\n");fclose(f);}
static double mapx(double x,double lo,double hi){return 110+740*(x-lo)/(hi-lo);}
static double mapy(double y,double lo,double hi){return 620-485*(y-lo)/(hi-lo);}
static void axes(FILE *f,double xmin,double xmax,double ymin,double ymax,const char *xl,const char *yl){
    if(xmax<=xmin)xmax=xmin+1;if(ymax<=ymin)ymax=ymin+1;
    for(int i=0;i<=5;i++){
        double x=110+740*i/5.0,y=135+485*i/5.0;
        fprintf(f,"<line x1=\"%.2f\" y1=\"135\" x2=\"%.2f\" y2=\"620\" stroke=\"#e4eaf1\"/>\n",x,x);
        fprintf(f,"<line x1=\"110\" y1=\"%.2f\" x2=\"850\" y2=\"%.2f\" stroke=\"#e4eaf1\"/>\n",y,y);
        if(strcmp(xl,"Steps per period")!=0)
            fprintf(f,"<text x=\"%.2f\" y=\"645\" font-size=\"13\" text-anchor=\"middle\">%.4g</text>\n",x,xmin+(xmax-xmin)*i/5);
        fprintf(f,"<text x=\"99\" y=\"%.2f\" font-size=\"13\" text-anchor=\"end\">%.4g</text>\n",y+4,ymax-(ymax-ymin)*i/5);
    }
    if(strcmp(xl,"Steps per period")==0){
        const int ticks[]={50,100,200,400,800,1600};
        for(int k=0;k<6;k++){
            double xx=mapx(ticks[k],xmin,xmax);
            fprintf(f,"<line x1=\"%.2f\" y1=\"135\" x2=\"%.2f\" y2=\"620\" stroke=\"#e4eaf1\"/>\n",xx,xx);
            fprintf(f,"<text x=\"%.2f\" y=\"645\" font-size=\"13\" text-anchor=\"middle\">%d</text>\n",xx,ticks[k]);
        }
    }
    fprintf(f,"<rect x=\"110\" y=\"135\" width=\"740\" height=\"485\" fill=\"none\" stroke=\"#475467\"/>\n");
    fprintf(f,"<text x=\"480\" y=\"690\" font-size=\"18\" text-anchor=\"middle\">%s</text>\n",xl);
    fprintf(f,"<text x=\"24\" y=\"375\" font-size=\"18\" text-anchor=\"middle\" transform=\"rotate(-90 24 375)\">%s</text>\n",yl);
}
static void polyline(FILE *f,const XY *a,int n,double xmin,double xmax,double ymin,double ymax,int color){
    fprintf(f,"<polyline fill=\"none\" stroke=\"%s\" stroke-width=\"2.7\" points=\"",palette[color%6]);
    for(int i=0;i<n;i++)if(isfinite(a[i].x)&&isfinite(a[i].y))
        fprintf(f,"%.2f,%.2f ",mapx(a[i].x,xmin,xmax),mapy(a[i].y,ymin,ymax));
    fprintf(f,"\"/>\n");
}
static void legend(FILE *f,int i,const char *name){
    int y=165+36*i;
    fprintf(f,"<line x1=\"880\" y1=\"%d\" x2=\"915\" y2=\"%d\" stroke=\"%s\" stroke-width=\"4\"/>\n",y,y,palette[i%6]);
    fprintf(f,"<text x=\"926\" y=\"%d\" font-size=\"15\">%s</text>\n",y+5,name);
}
static FILE *chart(const char *filename,const char *title,const char *subtitle,double xmin,double xmax,double ymin,double ymax,const char *xl,const char *yl){
    FILE *f=fopen(filename,"w");if(!f){perror(filename);return NULL;}
    svg_begin(f,title,subtitle);axes(f,xmin,xmax,ymin,ymax,xl,yl);return f;
}
static void trajectory_chart(const char *filename,const char *title,const char *subtitle,Result *res[],int n,const char *names[]){
    double xmin=0,xmax=0,ymin=0,ymax=0;
    for(int j=0;j<n;j++)for(int i=0;i<res[j]->n;i++){
        Point p=res[j]->p[i];if(p.x<xmin)xmin=p.x;if(p.x>xmax)xmax=p.x;
        if(p.y<ymin)ymin=p.y;if(p.y>ymax)ymax=p.y;
    }
    double dx=fmax(1,(xmax-xmin)*0.09),dy=fmax(1,(ymax-ymin)*0.12);
    xmin-=dx;xmax+=dx;ymin-=dy;ymax+=dy;
    FILE *f=chart(filename,title,subtitle,xmin,xmax,ymin,ymax,"x, mm","y, mm");if(!f)return;
    for(int j=0;j<n;j++){
        XY *p=malloc(res[j]->n*sizeof(*p));if(!p)break;
        for(int i=0;i<res[j]->n;i++)p[i]=(XY){res[j]->p[i].x,res[j]->p[i].y};
        polyline(f,p,res[j]->n,xmin,xmax,ymin,ymax,j);free(p);legend(f,j,names[j]);
    }
    svg_end(f);
}
static void line_chart(const char *filename,const char *title,const char *subtitle,const char *xl,const char *yl,const XY data[][64],const int sizes[],int n,const char *names[],int logy){
    double xmin=INFINITY,xmax=-INFINITY,ymin=INFINITY,ymax=-INFINITY;
    for(int j=0;j<n;j++)for(int i=0;i<sizes[j];i++){
        double x=data[j][i].x,y=data[j][i].y;
        if(logy)y=log10(fmax(y,1e-16));
        if(x<xmin)xmin=x;if(x>xmax)xmax=x;if(y<ymin)ymin=y;if(y>ymax)ymax=y;
    }
    double dx=fmax((xmax-xmin)*0.07,0.001),dy=fmax((ymax-ymin)*0.14,0.1);
    FILE *f=chart(filename,title,subtitle,xmin-dx,xmax+dx,ymin-dy,ymax+dy,xl,yl);if(!f)return;
    for(int j=0;j<n;j++){
        XY temp[64];for(int i=0;i<sizes[j];i++)temp[i]=(XY){data[j][i].x,logy?log10(fmax(data[j][i].y,1e-16)):data[j][i].y};
        polyline(f,temp,sizes[j],xmin-dx,xmax+dx,ymin-dy,ymax+dy,j);
        for(int i=0;i<sizes[j];i++)fprintf(f,"<circle cx=\"%.2f\" cy=\"%.2f\" r=\"4\" fill=\"%s\"/>\n",mapx(temp[i].x,xmin-dx,xmax+dx),mapy(temp[i].y,ymin-dy,ymax+dy),palette[j%6]);
        legend(f,j,names[j]);
    }
    svg_end(f);
}
static void detector_chart(const char *filename,const double hits[5]){
    double lo=INFINITY,hi=-INFINITY;
    for(int i=0;i<5;i++){if(hits[i]<lo)lo=hits[i];if(hits[i]>hi)hi=hits[i];}
    double margin=fmax(1.0,(hi-lo)*0.18);
    FILE *f=chart(filename,"Detector positions","Five distinct isotope hits; no connecting lines",
                  19.4,24.6,lo-margin,hi+margin,"Mass, u","Detector x, mm");
    if(!f)return;
    for(int i=0;i<5;i++){
        double x=mapx(20+i,19.4,24.6),y=mapy(hits[i],lo-margin,hi+margin);
        fprintf(f,"<circle cx=\"%.3f\" cy=\"%.3f\" r=\"7\" fill=\"%s\"/>\n",x,y,palette[i]);
        fprintf(f,"<text x=\"%.3f\" y=\"%.3f\" font-size=\"15\" text-anchor=\"middle\">%.4f</text>\n",x,y-13,hits[i]);
    }
    svg_end(f);
}
static void save_sweep(const char *name,const char *header,const XY d[][64],const int *sizes,int n){
    FILE *f=fopen(name,"w");if(!f)return;fprintf(f,"series,%s\n",header);
    for(int j=0;j<n;j++)for(int i=0;i<sizes[j];i++)fprintf(f,"%d,%.12g,%.12g\n",j,d[j][i].x,d[j][i].y);
    fclose(f);
}
static void file_name(char *dst,size_t cap,const char *dir,const char *base){snprintf(dst,cap,"%s/%s",dir,base);}
static int ensure_dir(const char *dir){
    char cmd[512];if(strlen(dir)>300||strchr(dir,'\''))return 0;
    snprintf(cmd,sizeof cmd,"mkdir -p '%s'",dir);return system(cmd)==0;
}
static void report_base(Params p,const char *dir){
    char path[512];file_name(path,sizeof path,dir,"summary.csv");FILE *f=fopen(path,"w");
    if(!f){perror(path);return;}
    fprintf(f,"mass_u,radius_mm,detector_mm,analytic_mm,error_mm,energy_relative\n");
    Result *r[5];const char *names[]={"20 u","21 u","22 u","23 u","24 u"};
    for(int i=0;i<5;i++){
        Params x=p;x.mass=20+i;char base[80];snprintf(base,sizeof base,"trajectory_%d.csv",20+i);
        file_name(path,sizeof path,dir,base);r[i]=simulate(x,RK4,2000,0.55,path);
        if(!r[i]){fprintf(stderr,"Simulation failed\n");exit(1);}
        double radius=x.mass*AMU*x.v/(Q*x.B)*1e3;
        double exacthit=2*radius*cos(x.angle*PI/180);
        fprintf(f,"%d,%.12g,%.12g,%.12g,%.12g,%.12g\n",20+i,radius,r[i]->hit,exacthit,fabs(r[i]->hit-exacthit),r[i]->energy);
        printf("m=%d u: r=%.4f mm; detector x=%.4f mm; exact=%.4f mm\n",20+i,radius,r[i]->hit,exacthit);
    }
    fclose(f);
    file_name(path,sizeof path,dir,"01_trajectories.svg");trajectory_chart(path,"Isotope trajectories","RK4; B along z; E = 0; masses 20-24 u",r,5,names);
    XY data[6][64]={0};int sizes[6]={0};
    for(int i=0;i<5;i++){data[0][i]=(XY){20+i,r[i]->hit};sizes[0]++;}
    double hits[5];for(int i=0;i<5;i++)hits[i]=r[i]->hit;
    file_name(path,sizeof path,dir,"02_detector.svg");detector_chart(path,hits);
    for(int i=0;i<5;i++)free(r[i]);
}
static void sweeps(Params p,const char *dir){
    char path[512];XY data[6][64]={0};int sizes[6]={0};
    const char *single[]={"m = 22 u"};
    for(int i=0;i<11;i++){
        Params x=p;x.mass=22;x.B=p.B*(0.5+0.1*i);
        Result *r=simulate(x,RK4,1000,0.75,NULL);data[0][i]=(XY){x.B,r->hit};free(r);sizes[0]++;
    }
    file_name(path,sizeof path,dir,"03_magnetic_field.svg");line_chart(path,"Magnetic field influence","Other parameters fixed","B, T","Detector x, mm",data,sizes,1,single,0);
    file_name(path,sizeof path,dir,"03_magnetic_field.csv");save_sweep(path,"B_T,detector_mm",data,sizes,1);
    memset(data,0,sizeof data);memset(sizes,0,sizeof sizes);
    for(int i=0;i<11;i++){
        Params x=p;x.mass=22;x.v=p.v*(0.5+0.1*i);
        Result *r=simulate(x,RK4,1000,0.75,NULL);data[0][i]=(XY){x.v,r->hit};free(r);sizes[0]++;
    }
    file_name(path,sizeof path,dir,"04_speed.svg");line_chart(path,"Initial speed influence","Other parameters fixed","v, m/s","Detector x, mm",data,sizes,1,single,0);
    file_name(path,sizeof path,dir,"04_speed.csv");save_sweep(path,"v_m_s,detector_mm",data,sizes,1);
    memset(data,0,sizeof data);memset(sizes,0,sizeof sizes);
    for(int i=0;i<13;i++){
        Params x=p;x.mass=22;x.angle=-30+5*i;
        Result *r=simulate(x,RK4,1000,0.75,NULL);data[0][i]=(XY){x.angle,r->hit};free(r);sizes[0]++;
    }
    file_name(path,sizeof path,dir,"05_angle.svg");line_chart(path,"Entry angle influence","Angle measured from +y toward +x","Angle, degrees","Detector x, mm",data,sizes,1,single,0);
    file_name(path,sizeof path,dir,"05_angle.csv");save_sweep(path,"angle_deg,detector_mm",data,sizes,1);
    Params e=p;e.mass=22;e.Ex=0;e.Ey=0;e.Ez=0;e.angle=0;
    Result *curves[5];const char *enames[]={"E = 0 V/m","Ex = 0.05 V/m","Ex = 0.10 V/m","Ex = 0.20 V/m","Ey = 0.10 V/m"};
    for(int i=0;i<5;i++){
        Params x=e;if(i==1)x.Ex=0.05;if(i==2)x.Ex=0.10;if(i==3)x.Ex=0.20;if(i==4)x.Ey=0.10;
        curves[i]=simulate(x,RK4,1600,1.5,NULL);
    }
    file_name(path,sizeof path,dir,"06_electric_field.svg");trajectory_chart(path,"Electric field influence","m = 22 u; B along z; Lorentz force",curves,5,enames);
    for(int i=0;i<5;i++)free(curves[i]);
    Result *zcurves[3];const char *znames[]={"vz = 0 m/s","vz = 250 m/s","vz = 500 m/s"};
    for(int i=0;i<3;i++){Params x=e;x.vz=250*i;zcurves[i]=simulate(x,RK4,1600,1.0,NULL);}
    memset(data,0,sizeof data);memset(sizes,0,sizeof sizes);
    for(int j=0;j<3;j++)for(int i=0;i<zcurves[j]->n;i+=32){
        Point t=zcurves[j]->p[i];data[j][sizes[j]++]=(XY){t.x,t.z};
    }
    file_name(path,sizeof path,dir,"07_longitudinal_speed.svg");line_chart(path,"Longitudinal velocity influence","Projection x-z of helical motion","x, mm","z, mm",data,sizes,3,znames,0);
    for(int i=0;i<3;i++)free(zcurves[i]);
    memset(data,0,sizeof data);memset(sizes,0,sizeof sizes);
    int ns[]={50,100,200,400,800,1600};const char *schemes[]={"RK4","Boris","Euler"};
    file_name(path,sizeof path,dir,"convergence.csv");FILE *f=fopen(path,"w");
    if(f)fprintf(f,"scheme,steps_per_period,max_error_mm,max_energy_relative\n");
    for(int j=0;j<3;j++)for(int i=0;i<6;i++){
        Result *r=simulate(e,j,ns[i],1.0,NULL);
        data[j][sizes[j]++]=(XY){ns[i],r->err};
        if(f)fprintf(f,"%s,%d,%.12g,%.12g\n",schemes[j],ns[i],r->err,r->energy);
        free(r);
    }
    if(f)fclose(f);
    file_name(path,sizeof path,dir,"08_convergence.svg");line_chart(path,"Convergence of numerical methods","Maximum position error over one period; log10(error in mm)","Steps per period","log10(error, mm)",data,sizes,3,schemes,1);
    memset(data,0,sizeof data);memset(sizes,0,sizeof sizes);
    for(int j=0;j<3;j++)for(int i=0;i<6;i++){
        Result *r=simulate(e,j,ns[i],1.0,NULL);data[j][sizes[j]++]=(XY){ns[i],r->energy};free(r);
    }
    file_name(path,sizeof path,dir,"09_energy.svg");line_chart(path,"Kinetic energy conservation","E = 0; maximum relative change; log10(error)","Steps per period","log10(relative energy change)",data,sizes,3,schemes,1);
    Result *compare[2];const char *cnames[]={"RK4 numerical","Analytic"};
    compare[0]=simulate(e,RK4,400,0.5,NULL);compare[1]=calloc(1,sizeof(Result));
    compare[1]->n=compare[0]->n;
    for(int i=0;i<compare[1]->n;i++){
        double t=compare[0]->p[i].t;State s=exact(e,t);
        compare[1]->p[i]=(Point){t,s.x*1e3,s.y*1e3,s.z*1e3};
    }
    file_name(path,sizeof path,dir,"10_analytic_vs_numeric.svg");trajectory_chart(path,"Analytic vs numerical","m = 22 u; RK4; E = 0",compare,2,cnames);
    free(compare[0]);free(compare[1]);
}
static void write_readme(const char *dir){
    char path[512];file_name(path,sizeof path,dir,"README.txt");FILE *f=fopen(path,"w");if(!f)return;
    fprintf(f,"Isotope separation, ITMO task 6.5\n\nRun on macOS:\ncd ~/Downloads\nclang -std=c11 -O2 main.c -o isotopes -lm\n./isotopes\n\nChange inputs:\n./isotopes --B 0.02 --v 1200 --angle 10\n./isotopes --Ex 0.1 --vz 200\n./isotopes --no-open\n\nAll output is in isotope_results/. The program generates 10 SVG plots, CSV data, and a summary.\nThe baseline plots show the specified mass separation. Parameter sweeps use m=22 u and fixed baseline values unless stated otherwise.\nCSV files contain raw calculated data; SVG files are directly viewable plots.\nMethods: Euler, Boris and RK4; analytic solution used only for E=0.\n");fclose(f);
}
int main(int argc,char **argv){
    Params p={22,0.01,0,0,0,1000,0,0};int open_graphs=1;
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"--no-open")){open_graphs=0;continue;}
        if(i+1>=argc){fprintf(stderr,"Missing value after %s\n",argv[i]);return 1;}
        char *end=NULL;double v=strtod(argv[i+1],&end);
        if(end==argv[i+1]||*end||!isfinite(v)){fprintf(stderr,"Invalid value: %s\n",argv[i+1]);return 1;}
        if(!strcmp(argv[i],"--B"))p.B=v;
        else if(!strcmp(argv[i],"--v"))p.v=v;
        else if(!strcmp(argv[i],"--angle"))p.angle=v;
        else if(!strcmp(argv[i],"--Ex"))p.Ex=v;
        else if(!strcmp(argv[i],"--Ey"))p.Ey=v;
        else if(!strcmp(argv[i],"--Ez"))p.Ez=v;
        else if(!strcmp(argv[i],"--vz"))p.vz=v;
        else{fprintf(stderr,"Unknown argument: %s\n",argv[i]);return 1;}
        i++;
    }
    if(p.B<=0||p.v<=0||fabs(p.angle)>=40||fabs(p.Ex)>2||fabs(p.Ey)>2||fabs(p.Ez)>2){
        fprintf(stderr,"Expected B>0, v>0, |angle|<40, |E|<=2 V/m\n");return 1;
    }
    const char *dir="isotope_results";if(!ensure_dir(dir)){fprintf(stderr,"Cannot create output folder\n");return 1;}
    printf("ITMO 6.5 | B=%.6g T, v=%.6g m/s, angle=%.3g deg\n",p.B,p.v,p.angle);
    if(!zero_e(p))printf("Note: baseline detector formulas assume E=0; electric field is studied separately.\n");
    Params base=p;base.Ex=base.Ey=base.Ez=0;base.vz=0;
    report_base(base,dir);sweeps(base,dir);write_readme(dir);
    if(!zero_e(p)||p.vz!=0){
        char path[512];file_name(path,sizeof path,dir,"custom_trajectory.csv");
        Result *custom=simulate(p,RK4,1600,1.0,path);
        if(custom){
            Result *arr[]={custom};const char *names[]={"Custom inputs"};
            file_name(path,sizeof path,dir,"11_custom_trajectory.svg");
            trajectory_chart(path,"Custom parameters","Includes requested E and vz",arr,1,names);
            if(p.vz!=0){
                XY data[6][64]={0};int sizes[6]={0};
                for(int i=0;i<custom->n;i+=32){Point t=custom->p[i];data[0][sizes[0]++]=(XY){t.x,t.z};}
                file_name(path,sizeof path,dir,"12_custom_xz.svg");
                line_chart(path,"Custom longitudinal motion","x-z projection","x, mm","z, mm",data,sizes,1,names,0);
            }
            free(custom);
        }
    }
    printf("Created 10 standard SVG graphs and CSV files in %s/\n",dir);
#ifdef __APPLE__
    if(open_graphs)system("open isotope_results/01_trajectories.svg isotope_results/02_detector.svg");
#else
    (void)open_graphs;
#endif
    return 0;
}
