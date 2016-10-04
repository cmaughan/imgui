#pragma once

#include <vector>

typedef unsigned int color32;

const color32 COLOR32_BLACK     = 0xff000000;   // Black 
const color32 COLOR32_WHITE     = 0xffffffff;   // White 
const color32 COLOR32_ZERO      = 0x00000000;   // Zero 
const color32 COLOR32_RED       = 0xffff0000;   // Red 
const color32 COLOR32_GREEN     = 0xff00ff00;   // Green 
const color32 COLOR32_BLUE      = 0xff0000ff;   // Blue 

template <typename _T> inline const _T& TClamp(const _T& _X, const _T& _Limit1, const _T& _Limit2)
{
    if( _Limit1<_Limit2 )
        return (_X<=_Limit1) ? _Limit1 : ( (_X>=_Limit2) ? _Limit2 : _X );
    else
        return (_X<=_Limit2) ? _Limit2 : ( (_X>=_Limit1) ? _Limit1 : _X );
}

inline color32 Color32FromARGBi(int _A, int _R, int _G, int _B)
{
    return (((color32)TClamp(_A, 0, 255))<<24) | (((color32)TClamp(_R, 0, 255))<<16) | (((color32)TClamp(_G, 0, 255))<<8) | ((color32)TClamp(_B, 0, 255));
}

inline color32 Color32FromARGBf(float _A, float _R, float _G, float _B)
{
    return (((color32)TClamp(_A*256.0f, 0.0f, 255.0f))<<24) | (((color32)TClamp(_R*256.0f, 0.0f, 255.0f))<<16) | (((color32)TClamp(_G*256.0f, 0.0f, 255.0f))<<8) | ((color32)TClamp(_B*256.0f, 0.0f, 255.0f));
}

inline void Color32ToARGBi(color32 _Color, int *_A, int *_R, int *_G, int *_B)
{
    if(_A) *_A = (_Color>>24)&0xff;
    if(_R) *_R = (_Color>>16)&0xff;
    if(_G) *_G = (_Color>>8)&0xff;
    if(_B) *_B = _Color&0xff;
}

inline void Color32ToARGBf(color32 _Color, float *_A, float *_R, float *_G, float *_B)
{
    if(_A) *_A = (1.0f/255.0f)*float((_Color>>24)&0xff);
    if(_R) *_R = (1.0f/255.0f)*float((_Color>>16)&0xff);
    if(_G) *_G = (1.0f/255.0f)*float((_Color>>8)&0xff);
    if(_B) *_B = (1.0f/255.0f)*float(_Color&0xff);
}

inline void ColorRGBToHLSf(float _R, float _G, float _B, float *_Hue, float *_Light, float *_Saturation)
{
    // Compute HLS from RGB. The r,g,b triplet is between [0,1], 
    // hue is between [0,360], light and saturation are [0,1].

    float rnorm, gnorm, bnorm, minval, maxval, msum, mdiff, r, g, b;
    r = g = b = 0;
    if(_R>0) r = _R; if(r>1) r = 1;
    if(_G>0) g = _G; if(g>1) g = 1;
    if(_B>0) b = _B; if(b>1) b = 1;

    minval = r;
    if(g<minval) minval = g;
    if(b<minval) minval = b;
    maxval = r;
    if(g>maxval) maxval = g;
    if(b>maxval) maxval = b;

    rnorm = gnorm = bnorm = 0;
    mdiff = maxval - minval;
    msum  = maxval + minval;
    float l = 0.5f * msum;
    if(_Light) 
        *_Light = l;
    if(maxval!=minval) 
    {
        rnorm = (maxval - r)/mdiff;
        gnorm = (maxval - g)/mdiff;
        bnorm = (maxval - b)/mdiff;
    } 
    else 
    {
        if(_Saturation)
            *_Saturation = 0;
        if(_Hue)
            *_Hue = 0;
        return;
    }

    if(_Saturation)
    {
        if(l<0.5f)
            *_Saturation = mdiff/msum;
        else
            *_Saturation = mdiff/(2.0f - msum);
    }

    if(_Hue)
    {
        if(r==maxval)
            *_Hue = 60.0f * (6.0f + bnorm - gnorm);
        else if(g==maxval)
            *_Hue = 60.0f * (2.0f + rnorm - bnorm);
        else
            *_Hue = 60.0f * (4.0f + gnorm - rnorm);

        if(*_Hue>360.0f)
            *_Hue -= 360.0f;
    }
}

inline void ColorRGBToHLSi(int _R, int _G, int _B, int *_Hue, int *_Light, int *_Saturation)
{
    float h, l, s;
    ColorRGBToHLSf((1.0f/255.0f)*float(_R), (1.0f/255.0f)*float(_G), (1.0f/255.0f)*float(_B), &h, &l, &s);
    if(_Hue)        *_Hue       = (int)TClamp(h*(256.0f/360.0f), 0.0f, 255.0f);
    if(_Light)      *_Light     = (int)TClamp(l*256.0f, 0.0f, 255.0f);
    if(_Saturation) *_Saturation= (int)TClamp(s*256.0f, 0.0f, 255.0f);
}

inline void ColorHLSToRGBf(float _Hue, float _Light, float _Saturation, float *_R, float *_G, float *_B)
{
    // Compute RGB from HLS. The light and saturation are between [0,1]
    // and hue is between [0,360]. The returned r,g,b triplet is between [0,1].

    // a local auxiliary function
    struct CLocal
    {
        static float HLSToRGB(float _Rn1, float _Rn2, float _Huei)
        {
            float hue = _Huei;
            if(hue>360) hue = hue - 360;
            if(hue<0)   hue = hue + 360;
            if(hue<60 ) return _Rn1 + (_Rn2-_Rn1)*hue/60;
            if(hue<180) return _Rn2;
            if(hue<240) return _Rn1 + (_Rn2-_Rn1)*(240-hue)/60;
            return _Rn1;
        }
    };

    float rh, rl, rs, rm1, rm2;
    rh = rl = rs = 0;
    if(_Hue>0)        rh = _Hue;        if(rh>360) rh = 360;
    if(_Light>0)      rl = _Light;      if(rl>1)   rl = 1;
    if(_Saturation>0) rs = _Saturation; if(rs>1)   rs = 1;

    if(rl<=0.5f)
        rm2 = rl*(1.0f + rs);
    else
        rm2 = rl + rs - rl*rs;
    rm1 = 2.0f*rl - rm2;

    if(!rs) 
    { 
        if(_R) *_R = rl; 
        if(_G) *_G = rl; 
        if(_B) *_B = rl; 
    }
    else
    {
        if(_R) *_R = CLocal::HLSToRGB(rm1, rm2, rh+120);
        if(_G) *_G = CLocal::HLSToRGB(rm1, rm2, rh);
        if(_B) *_B = CLocal::HLSToRGB(rm1, rm2, rh-120);
    }
}

inline void ColorHLSToRGBi(int _Hue, int _Light, int _Saturation, int *_R, int *_G, int *_B)
{
    float r, g, b;
    ColorHLSToRGBf((360.0f/255.0f)*float(_Hue), (1.0f/255.0f)*float(_Light), (1.0f/255.0f)*float(_Saturation), &r, &g, &b);
    if(_R) *_R = (int)TClamp(r*256.0f, 0.0f, 255.0f);
    if(_G) *_G = (int)TClamp(g*256.0f, 0.0f, 255.0f);
    if(_B) *_B = (int)TClamp(b*256.0f, 0.0f, 255.0f);
}


inline color32 ColorBlend(color32 _Color1, color32 _Color2, float _S)
{
    float a1, r1, g1, b1, a2, r2, g2, b2;
    Color32ToARGBf(_Color1, &a1, &r1, &g1, &b1);
    Color32ToARGBf(_Color2, &a2, &r2, &g2, &b2);
    float t = 1.0f-_S;
    return Color32FromARGBf(t*a1+_S*a2, t*r1+_S*r2, t*g1+_S*g2, t*b1+_S*b2);
}

struct IMGUI_API ImQuat
{
    ImQuat();
    double               Qx, Qy, Qz, Qs;    // Quat value
    double               Vx, Vy, Vz, Angle; // Not used
    double               Dx, Dy, Dz;        // Dir value set when used as a direction
    bool                 m_AAMode;          // Axis & angle mode -> disabled
    bool                 m_ShowVal;         // Display values
    bool                 m_IsFloat;         // Quat/Dir uses floats
    bool                 m_IsDir;           // Mapped to a dir vector instead of a quat
    double               m_Dir[3];          // If not zero, display one direction vector
    color32              m_DirColor;        // Direction vector color
    float                m_Permute[3][3];   // Permute frame axis
    static std::vector<float>   s_SphTri;
    static std::vector<color32> s_SphCol;
    static std::vector<int>     s_SphTriProj;
    static std::vector<color32> s_SphColLight;
    static std::vector<float>   s_ArrowTri[4];
    static std::vector<int>     s_ArrowTriProj[4];
    static std::vector<float>   s_ArrowNorm[4];
    static std::vector<color32> s_ArrowColLight[4];
    
    enum EArrowParts     { ARROW_CONE, ARROW_CONE_CAP, ARROW_CYL, ARROW_CYL_CAP };
    static void          CreateSphere();
    static void          CreateArrow();
    void                 ConvertToAxisAngle();
    void                 ConvertFromAxisAngle();
    static void          ApplyQuat(float *outX, float *outY, float *outZ, float x, float y, float z, float qx, float qy, float qz, float qs);
    static void          QuatFromDir(double *outQx, double *outQy, double *outQz, double *outQs, double dx, double dy, double dz);
    inline void          Permute(float *outX, float *outY, float *outZ, float x, float y, float z);
    inline void          PermuteInv(float *outX, float *outY, float *outZ, float x, float y, float z);
    inline void          Permute(double *outX, double *outY, double *outZ, double x, double y, double z);
    inline void          PermuteInv(double *outX, double *outY, double *outZ, double x, double y, double z);
    /*
    static void IMGUI_API CopyVarFromExtCB(void *_VarValue, const void *_ExtValue, unsigned int _ExtMemberIndex, void *_ClientData);
    static void IMGUI_API CopyVarToExtCB(const void *_VarValue, void *_ExtValue, unsigned int _ExtMemberIndex, void *_ClientData);
    static void IMGUI_API SummaryCB(char *_SummaryString, size_t _SummaryMaxLength, const void *_ExtValue, void *_ClientData);
    static void IMGUI_API DrawCB(int _W, int _H, void *_ExtValue, void *_ClientData);
    static bool IMGUI_API MouseMotionCB(int _MouseX, int _MouseY, int _W, int _H, void *_StructExtValue, void *_ClientData);
    static bool IMGUI_API MouseButtonCB(TwMouseButtonID _Button, bool _Pressed, int _MouseX, int _MouseY, int _W, int _H, void *_StructExtValue, void *_ClientData);
    static void IMGUI_API MouseLeaveCB(void *_StructExtValue, void *_ClientData);
*/    
    bool                 m_Highlighted;
    bool                 m_Rotating;
    double               m_OrigQuat[4];
    float                m_OrigX, m_OrigY;
    double               m_PrevX, m_PrevY;

};

bool IMGUI_API ImGui_Orient(char* label, ImQuat& quat);

static inline void QuatMult(double *out, const double *q1, const double *q2)
{
    out[0] = q1[3]*q2[0] + q1[0]*q2[3] + q1[1]*q2[2] - q1[2]*q2[1];
    out[1] = q1[3]*q2[1] + q1[1]*q2[3] + q1[2]*q2[0] - q1[0]*q2[2];
    out[2] = q1[3]*q2[2] + q1[2]*q2[3] + q1[0]*q2[1] - q1[1]*q2[0];
    out[3] = q1[3]*q2[3] - (q1[0]*q2[0] + q1[1]*q2[1] + q1[2]*q2[2]);
}

static inline void QuatFromAxisAngle(double *out, const double *axis, double angle)
{
    double n = axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2];
    if( fabs(n)>DBL_EPSILON )
    {
        double f = 0.5*angle;
        out[3] = cos(f);
        f = sin(f)/sqrt(n);
        out[0] = axis[0]*f;
        out[1] = axis[1]*f;
        out[2] = axis[2]*f;
    }
    else
    {
        out[3] = 1.0;
        out[0] = out[1] = out[2] = 0.0;
    }
}

static inline void Vec3Cross(double *out, const double *a, const double *b)
{
    out[0] = a[1]*b[2]-a[2]*b[1];
    out[1] = a[2]*b[0]-a[0]*b[2];
    out[2] = a[0]*b[1]-a[1]*b[0];
}

static inline double Vec3Dot(const double *a, const double *b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline void Vec3RotY(float *x, float *y, float *z)
{
    (void)y;
    float tmp = *x;
    *x = - *z;
    *z = tmp;
}

static inline void Vec3RotZ(float *x, float *y, float *z)
{
    (void)z;
    float tmp = *x;
    *x = - *y;
    *y = tmp;
}

static inline float QuatD(int w, int h)
{
    return (float)std::min(abs(w), abs(h)) - 4;
}

static inline int QuatPX(float x, int w, int h)
{
    return (int)(x*0.5f*QuatD(w, h) + (float)w*0.5f + 0.5f);
}

static inline int QuatPY(float y, int w, int h)
{
    return (int)(-y*0.5f*QuatD(w, h) + (float)h*0.5f - 0.5f);
}

static inline float QuatIX(int x, int w, int h)
{
    return (2.0f*(float)x - (float)w - 1.0f)/QuatD(w, h);
}

static inline float QuatIY(int y, int w, int h)
{
    return (-2.0f*(float)y + (float)h - 1.0f)/QuatD(w, h);
}

enum Cull { CULL_NONE, CULL_CW, CULL_CCW };
const float  FLOAT_EPS     = 1.0e-7f;
const float  FLOAT_EPS_SQ  = 1.0e-14f;
const float  FLOAT_PI      = 3.14159265358979323846f;
const double DOUBLE_EPS    = 1.0e-14;
const double DOUBLE_EPS_SQ = 1.0e-28;
const double DOUBLE_PI     = 3.14159265358979323846;

inline double DegToRad(double degree) { return degree * (DOUBLE_PI/180.0); }
inline double RadToDeg(double radian) { return radian * (180.0/DOUBLE_PI); }
