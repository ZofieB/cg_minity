//------------------------------------------------------------------
#version 420 core
//------------------------------------------------------------------
// Ray tracer ver: 1.000
//------------------------------------------------------------------
in smooth vec3      ray_pos;    // ray start position
in smooth vec3      ray_dir;    // ray start direction
uniform float       n0;         // refractive index of camera origin
uniform int         fac_siz;    // square texture x,y resolution size
uniform int         fac_num;    // number of valid floats in texture
uniform sampler2D   fac_txr;    // scene mesh data texture
out layout(location=0) vec4 frag_col;
//---------------------------------------------------------------------------
//#define _debug_print
#define _reflect
#define _refract
//---------------------------------------------------------------------------
#ifdef _debug_print
in vec2 txt_pos;                // frag screen position <-1,+1>
uniform sampler2D txr_font;     // ASCII 32x8 characters font texture unit
uniform float txt_fxs,txt_fys;  // font/screen resolution ratio
const int _txtsiz=64;           // text buffer size
int txt[_txtsiz],txtsiz;        // text buffer and its actual size
vec4 txt_col=vec4(0.0,0.0,0.0,1.0); // color interface for txt_print()
bool _txt_col=false;            // is txt_col active?
void txt_decimal(vec2 v);       // print vec3 into txt
void txt_decimal(vec3 v);       // print vec3 into txt
void txt_decimal(vec4 v);       // print vec3 into txt
void txt_decimal(float x);      // print float x into txt
void txt_decimal(int x);        // print int x into txt
void txt_print(float x0,float y0);  // print txt at x0,y0 [chars]
#endif
//---------------------------------------------------------------------------
void main(void)
    {
    const vec3  light_dir=normalize(vec3(0.1,0.1,1.0));
    const float light_iamb=0.1;                 // dot offset
    const float light_idir=0.5;                 // directional light amplitude
    const vec3 back_col=vec3(0.2,0.2,0.2);      // background color

    const float _zero=1e-6;     // to avoid intrsection with start point of ray
    const int _fac_triangles=0; // r,g,b, refl,refr,n, type, triangle count, { x0,y0,z0,x1,y1,z1,x2,y2,z2 }
    const int _fac_spheres  =1; // r,g,b, refl,refr,n, type, sphere count,   { x,y,z,r }
    // ray scene intersection
    struct _ray
        {
        vec3 pos,dir,nor;
        vec3 col;
        float refl,refr;// reflection,refraction intensity coeficients
        float n0,n1,l;  // refaction index (start,end) , ray length
        int lvl,i0,i1;  // recursion level, reflect, refract
        };
    const int _lvls=5;
    const int _rays=(1<<_lvls)-1;
    _ray ray[_rays]; int rays;

    vec3 v0,v1,v2,pos;
    vec3 c,col;
    float refr,refl;
    float tt,t,n1,a;
    int i0,ii,num,id;

    // fac texture access
    vec2 st; int i,j; float ds=1.0/float(fac_siz-1);
    #define fac_get texture(fac_txr,st).r; st.s+=ds; i++; j++; if (j==fac_siz) { j=0; st.s=0.0; st.t+=ds; }
    // enque start ray
    ray[0].pos=ray_pos;
    ray[0].dir=normalize(ray_dir);
    ray[0].nor=vec3(0.0,0.0,0.0);
    ray[0].refl=0.0;
    ray[0].refr=0.0;
    ray[0].n0=n0;
    ray[0].n1=1.0;
    ray[0].l =0.0;
    ray[0].lvl=0;
    ray[0].i0=-1;
    ray[0].i1=-1;
    rays=1;

    // debug print area
    #ifdef _debug_print
    bool _dbg=false;
    float dbg_x0=45.0;
    float dbg_y0= 1.0;
    float dbg_xs=12.0;
    float dbg_ys=_rays+1.0;

    dbg_xs=40.0;
    dbg_ys=10;

    float x=0.5*(1.0+txt_pos.x)/txt_fxs; x-=dbg_x0;
    float y=0.5*(1.0-txt_pos.y)/txt_fys; y-=dbg_y0;
    // inside bbox?
    if ((x>=0.0)&&(x<=dbg_xs)
      &&(y>=0.0)&&(y<=dbg_ys))
        {
        // prints on
        _dbg=true;
        // preset debug ray
        ray[0].pos=vec3(0.0,0.0,0.0)*2.5;
        ray[0].dir=vec3(0.0,0.0,1.0);
        }
    #endif

    // loop all enqued rays
    for (i0=0;i0<rays;i0++)
        {
        // loop through all objects
        // find closest forward intersection between them and ray[i0]
        // strore it to ray[i0].(nor,col)
        // strore it to pos,n1
        t=tt=-1.0; ii=1; ray[i0].l=0.0;
        ray[i0].col=back_col;
        pos=ray[i0].pos; n1=n0;
        for (st=vec2(0.0,0.0),i=j=0;i<fac_num;)
            {
            c.r=fac_get;            // RGBA
            c.g=fac_get;
            c.b=fac_get;
            refl=fac_get;
            refr=fac_get;
            n1=fac_get;             // refraction index
            a=fac_get; id=int(a);   // object type
            a=fac_get; num=int(a);  // face count

            if (id==_fac_triangles)
             for (;num>0;num--)
                {
                v0.x=fac_get; v0.y=fac_get; v0.z=fac_get;
                v1.x=fac_get; v1.y=fac_get; v1.z=fac_get;
                v2.x=fac_get; v2.y=fac_get; v2.z=fac_get;
                vec3 e1,e2,n,p,q,r;
                float t,u,v,det,idet;
                //compute ray triangle intersection
                e1=v1-v0;
                e2=v2-v0;
                // Calculate planes normal vector
                p=cross(ray[i0].dir,e2);
                det=dot(e1,p);
                // Ray is parallel to plane
                if (abs(det)<1e-8) continue;
                idet=1.0/det;
                r=ray[i0].pos-v0;
                u=dot(r,p)*idet;
                if ((u<0.0)||(u>1.0)) continue;
                q=cross(r,e1);
                v=dot(ray[i0].dir,q)*idet;
                if ((v<0.0)||(u+v>1.0)) continue;
                t=dot(e2,q)*idet;
                if ((t>_zero)&&((t<=tt)||(ii!=0)))
                    {
                    ii=0; tt=t;
                    // store color,n ...
                    ray[i0].col=c;
                    ray[i0].refl=refl;
                    ray[i0].refr=refr;
                    // barycentric interpolate position
                    t=1.0-u-v;
                    pos=(v0*t)+(v1*u)+(v2*v);
                    // compute normal (store as dir for now)
                    e1=v1-v0;
                    e2=v2-v1;
                    ray[i0].nor=cross(e1,e2);
                    }
                }

            if (id==_fac_spheres)
             for (;num>0;num--)
                {
                float r;
                v0.x=fac_get; v0.y=fac_get; v0.z=fac_get; r=fac_get;
                // compute l0 length of ray(p0,dp) to intersection with sphere(v0,r)
                // where rr= r^-2
                float aa,bb,cc,dd,l0,l1,rr;
                vec3 p0,dp;
                p0=ray[i0].pos-v0;  // set sphere center to (0,0,0)
                dp=ray[i0].dir;
                rr = 1.0/(r*r);
                aa=2.0*rr*dot(dp,dp);
                bb=2.0*rr*dot(p0,dp);
                cc=    rr*dot(p0,p0)-1.0;
                dd=((bb*bb)-(2.0*aa*cc));
                if (dd<0.0) continue;
                dd=sqrt(dd);
                l0=(-bb+dd)/aa;
                l1=(-bb-dd)/aa;
                if (l0<0.0) l0=l1;
                if (l1<0.0) l1=l0;
                t=min(l0,l1); if (t<=_zero) t=max(l0,l1);
                if ((t>_zero)&&((t<=tt)||(ii!=0)))
                    {
                    ii=0; tt=t;
                    // store color,n ...
                    ray[i0].col=c;
                    ray[i0].refl=refl;
                    ray[i0].refr=refr;
                    // position,normal
                    pos=ray[i0].pos+(ray[i0].dir*t);
                    ray[i0].nor=pos-v0;
                    }
                }
            }
        ray[i0].l=tt;
        ray[i0].nor=normalize(ray[i0].nor);
        // split ray from pos and ray[i0].nor
        if ((ii==0)&&(ray[i0].lvl<_lvls-1))
            {
            t=dot(ray[i0].dir,ray[i0].nor);

            // reflect
            #ifdef _reflect
            if ((ray[i0].refl>_zero)&&(t<_zero))    // do not reflect inside objects
                {
                ray[i0].i0=rays;
                ray[rays]=ray[i0];
                ray[rays].lvl++;
                ray[rays].i0=-1;
                ray[rays].i1=-1;
                ray[rays].pos=pos;
                ray[rays].dir=ray[rays].dir-(2.0*t*ray[rays].nor);
                ray[rays].n0=ray[i0].n0;
                ray[rays].n1=ray[i0].n0;
                rays++;
                }
            #endif

            // refract
            #ifdef _refract
            if (ray[i0].refr>_zero)
                {
                ray[i0].i1=rays;
                ray[rays]=ray[i0];
                ray[rays].lvl++;
                ray[rays].i0=-1;
                ray[rays].i1=-1;
                ray[rays].pos=pos;

                t=dot(ray[i0].dir,ray[i0].nor);
                if (t>0.0)  // exit object
                    {
                    ray[rays].n0=ray[i0].n0;
                    ray[rays].n1=n0;
                    v0=-ray[i0].nor; t=-t;
                    }
                else{       // enter object
                    ray[rays].n0=n1;
                    ray[rays].n1=ray[i0].n0;
                    ray[i0  ].n1=n1;
                    v0=ray[i0].nor;
                    }
                n1=ray[i0].n0/ray[i0].n1;
                tt=1.0-(n1*n1*(1.0-t*t));
                if (tt>=0.0)
                    {
                    ray[rays].dir=(ray[i0].dir*n1)-(v0*((n1*t)+sqrt(tt)));
                    rays++;
                    }
                }
            #endif
            }
        else if (i0>0) // ignore last ray if nothing hit
            {
            ray[i0]=ray[rays-1];
            rays--; i0--;
            }
        }
    // back track ray intersections and compute output color col
    // lvl is sorted ascending so backtrack from end
    for (i0=rays-1;i0>=0;i0--)
        {
        // directional + ambient light
        t=abs(dot(ray[i0].nor,light_dir)*light_idir)+light_iamb;
        t*=1.0-ray[i0].refl-ray[i0].refr;
        ray[i0].col.rgb*=t;
        // reflect
        ii=ray[i0].i0;
        if (ii>=0) ray[i0].col.rgb+=ray[ii].col.rgb*ray[i0].refl;
        // refract
        ii=ray[i0].i1;
        if (ii>=0) ray[i0].col.rgb+=ray[ii].col.rgb*ray[i0].refr;
        }

    col=ray[0].col;

    // debug prints
    #ifdef _debug_print
/*
    if (_dbg)
        {
        txtsiz=0;
        txt_decimal(_lvls);
        txt[txtsiz]=' '; txtsiz++;
        txt_decimal(rays);
        txt[txtsiz]=' '; txtsiz++;
        txt_decimal(_rays);
        txt_print(dbg_x0,dbg_y0);

        for (ii=0;ii<rays;ii++)
            {
            txtsiz=0;
            txt_decimal(ray[ii].lvl);
            txt_print(dbg_x0,dbg_y0+ii+1);
            }

        for (ii=0,st=vec2(0.0,0.0),i=j=0;i<fac_num;ii++)
            {
            c.r=fac_get;            // RGBA
            txtsiz=0;
            txt_decimal(c.r);
            txt_print(dbg_x0,dbg_y0+ii+1);
            }
        if (_txt_col) col=txt_col.rgb;
        }
*/
    if (_dbg)
        {
        float x=dbg_x0,y=dbg_y0;
        vec3 a=vec3(1.0,2.0,3.0);
        vec3 b=vec3(5.0,6.0,7.0);
        txtsiz=0; txt_decimal(dot(a,b)); txt_print(x,y); y++;
        txtsiz=0; txt_decimal(cross(a,b)); txt_print(x,y); y++;
        if (_txt_col) col=txt_col.rgb;
        }
    #endif

    frag_col=vec4(col,1.0);
    }
//---------------------------------------------------------------------------
#ifdef _debug_print
//---------------------------------------------------------------------------
void txt_decimal(vec2 v)        // print vec2 into txt
    {
                      txt[txtsiz]='('; txtsiz++;
    txt_decimal(v.x); txt[txtsiz]=','; txtsiz++;
    txt_decimal(v.y); txt[txtsiz]=')'; txtsiz++;
    txt[txtsiz]=0;  // string terminator
    }
//---------------------------------------------------------------------------
void txt_decimal(vec3 v)        // print vec3 into txt
    {
                      txt[txtsiz]='('; txtsiz++;
    txt_decimal(v.x); txt[txtsiz]=','; txtsiz++;
    txt_decimal(v.y); txt[txtsiz]=','; txtsiz++;
    txt_decimal(v.z); txt[txtsiz]=')'; txtsiz++;
    txt[txtsiz]=0;  // string terminator
    }
//---------------------------------------------------------------------------
void txt_decimal(vec4 v)        // print vec4 into txt
    {
                      txt[txtsiz]='('; txtsiz++;
    txt_decimal(v.x); txt[txtsiz]=','; txtsiz++;
    txt_decimal(v.y); txt[txtsiz]=','; txtsiz++;
    txt_decimal(v.z); txt[txtsiz]=','; txtsiz++;
    txt_decimal(v.w); txt[txtsiz]=')'; txtsiz++;
    txt[txtsiz]=0;  // string terminator
    }
//---------------------------------------------------------------------------
void txt_decimal(float x)       // print float x into txt
    {
    int i,j,c;                  // l is size of string
    float y,a;
    const float base=10;
    // handle sign
    if (x<0.0) { txt[txtsiz]='-'; txtsiz++; x=-x; }
     else      { txt[txtsiz]='+'; txtsiz++; }
    // divide to int(x).fract(y) parts of number
    y=x; x=floor(x); y-=x;
    // handle integer part
    i=txtsiz;                   // start of integer part
    for (;txtsiz<_txtsiz;)
        {
        a=x;
        x=floor(x/base);
        a-=base*x;
        txt[txtsiz]=int(a)+'0'; txtsiz++;
        if (x<=0.0) break;
        }
    j=txtsiz-1;                 // end of integer part
    for (;i<j;i++,j--)          // reverse integer digits
        {
        c=txt[i]; txt[i]=txt[j]; txt[j]=c;
        }
    // handle fractional part
    for (txt[txtsiz]='.',txtsiz++;txtsiz<_txtsiz;)
        {
        y*=base;
        a=floor(y);
        y-=a;
        txt[txtsiz]=int(a)+'0'; txtsiz++;
        if (y<=0.0) break;
        }
    txt[txtsiz]=0;  // string terminator
    }
//---------------------------------------------------------------------------
void txt_decimal(int x)     // print int x into txt
    {
    int a,i,j,c;            // l is size of string
    const int base=10;
    // handle sign
    if (x<0.0) { txt[txtsiz]='-'; txtsiz++; x=-x; }
     else      { txt[txtsiz]='+'; txtsiz++; }
    // handle integer part
    i=txtsiz;               // start of integer part
    for (;txtsiz<_txtsiz;)
        {
        a=x;
        x/=base;
        a-=base*x;
        txt[txtsiz]=int(a)+'0'; txtsiz++;
        if (x<=0) break;
        }
    j=txtsiz-1;             // end of integer part
    for (;i<j;i++,j--)      // reverse integer digits
        {
        c=txt[i]; txt[i]=txt[j]; txt[j]=c;
        }
    txt[txtsiz]=0;  // string terminator
    }
//---------------------------------------------------------------------------
void txt_print(float x0,float y0)   // print txt at x0,y0 [chars]
    {
    int i;
    float x,y;
    // fragment position [chars] relative to x0,y0
    x=0.5*(1.0+txt_pos.x)/txt_fxs; x-=x0;
    y=0.5*(1.0-txt_pos.y)/txt_fys; y-=y0;
    // inside bbox?
    if ((x<0.0)||(x>float(txtsiz))||(y<0.0)||(y>1.0)) return;
    // get font texture position for target ASCII
    i=int(x);               // char index in txt
    x-=float(i);
    i=txt[i];
    x+=float(int(i&31));
    y+=float(int(i>>5));
    x/=32.0; y/=8.0;    // offset in char texture
    txt_col=texture(txr_font,vec2(x,y));
    _txt_col=true;
    }
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------